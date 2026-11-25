#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from xatlas.frontend_bundle import default_panel_specs, export_frontend_bundle, panel_strata


def run_cmd(args: list[str], *, env: dict[str, str]) -> None:
    print("[run] " + " ".join(args))
    subprocess.run(args, check=True, cwd=REPO_ROOT, env=env)


def repo_env() -> dict[str, str]:
    env = os.environ.copy()
    pythonpath = str(REPO_ROOT / "src")
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = pythonpath if not existing else pythonpath + os.pathsep + existing
    return env


def ensure_shared_prereqs(
    args: argparse.Namespace,
    *,
    env: dict[str, str],
    shared_dir: Path,
) -> tuple[Path, Path]:
    genes_out = shared_dir / "ensembl_x_genes.tsv"
    eqtl_index_out = shared_dir / "eqtl_chrX_index.tsv"

    if args.force_shared or not genes_out.exists():
        run_cmd(
            [
                sys.executable,
                "-B",
                str(REPO_ROOT / "scripts" / "03_build_x_gene_catalog.py"),
                "--out",
                str(genes_out),
            ],
            env=env,
        )
    else:
        print(f"[skip] shared genes exists: {genes_out}")

    if args.force_shared or not eqtl_index_out.exists():
        run_cmd(
            [
                sys.executable,
                "-B",
                str(REPO_ROOT / "scripts" / "04_prepare_eqtl_index.py"),
                "--outdir",
                str(shared_dir / "eqtl_catalogue"),
                "--index-out",
                str(eqtl_index_out),
            ],
            env=env,
        )
    else:
        print(f"[skip] shared eQTL index exists: {eqtl_index_out}")

    return genes_out, eqtl_index_out


def selected_traits_count(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        df = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return 0
    return int(len(df))


def stratum_workdirs(work_root: Path, panel_id: str, stratum: str) -> dict[str, Path]:
    root = work_root / panel_id / stratum
    return {
        "root": root,
        "interim": root / "interim",
        "processed": root / "processed",
        "release": root / "release",
        "extracts": root / "interim" / "panukb_x",
    }


def existing_panel_specs(release_root: Path) -> list:
    specs = []
    for spec in default_panel_specs(release_root):
        strata = panel_strata(spec)
        has_release = bool(strata)
        if not has_release and spec.release_dir is not None:
            has_release = (Path(spec.release_dir) / "trait_scores.tsv").exists()
        if has_release:
            specs.append(spec)
    return specs


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a stratified backend build for a single panel.")
    parser.add_argument("--panel-id", required=True)
    parser.add_argument("--panel-cfg", required=True)
    parser.add_argument("--work-root", default=os.path.expanduser("~/tmp/xchratlas_stratified"))
    parser.add_argument("--release-root", default=None)
    parser.add_argument("--frontend-outdir", default=None)
    parser.add_argument("--strata", nargs="+", default=["both_sexes", "females", "males"])
    parser.add_argument("--strict-pheno-sex", action="store_true")
    parser.add_argument("--allow-fallback-both-sexes", action="store_true")
    parser.add_argument("--min-curated-traits-per-stratum", type=int, default=1)
    parser.add_argument("--force-shared", action="store_true")
    parser.add_argument("--skip-frontend-bundle", action="store_true")
    args = parser.parse_args()

    env = repo_env()
    work_root = Path(args.work_root).expanduser()
    release_root = Path(args.release_root).expanduser() if args.release_root else work_root / "release"
    frontend_outdir = Path(args.frontend_outdir).expanduser() if args.frontend_outdir else release_root / "frontend"
    shared_dir = work_root / "shared"
    shared_dir.mkdir(parents=True, exist_ok=True)

    genes_out, eqtl_index_out = ensure_shared_prereqs(args, env=env, shared_dir=shared_dir)

    built_strata: list[str] = []
    for stratum in args.strata:
        dirs = stratum_workdirs(work_root, args.panel_id, stratum)
        dirs["interim"].mkdir(parents=True, exist_ok=True)
        dirs["processed"].mkdir(parents=True, exist_ok=True)
        release_outdir = release_root / "stratified" / args.panel_id / stratum

        select_cmd = [
            sys.executable,
            "-B",
            str(REPO_ROOT / "scripts" / "01_select_seed_traits.py"),
            "--panel",
            args.panel_cfg,
            "--outdir",
            str(dirs["interim"]),
            "--pheno-sex",
            stratum,
            "--selection-audit-out",
            str(dirs["interim"] / "selection_audit.tsv"),
        ]
        if args.strict_pheno_sex:
            select_cmd.append("--strict-pheno-sex")
        if args.allow_fallback_both_sexes:
            select_cmd.append("--allow-fallback-both-sexes")
        run_cmd(select_cmd, env=env)

        selected_path = dirs["interim"] / "selected_traits.tsv"
        n_selected = selected_traits_count(selected_path)
        if n_selected == 0:
            print(f"[skip] {args.panel_id}/{stratum}: no selected traits")
            continue

        run_cmd(
            [
                sys.executable,
                "-B",
                str(REPO_ROOT / "scripts" / "02_extract_panukb_x.py"),
                "--selected",
                str(selected_path),
                "--outdir",
                str(dirs["extracts"]),
            ],
            env=env,
        )
        run_cmd(
            [
                sys.executable,
                "-B",
                str(REPO_ROOT / "scripts" / "05_call_x_loci.py"),
                "--extract-dir",
                str(dirs["extracts"]),
                "--selected",
                str(selected_path),
                "--out-loci",
                str(dirs["processed"] / "x_loci.tsv.gz"),
                "--out-trait-summary",
                str(dirs["processed"] / "x_trait_locus_summary.tsv"),
            ],
            env=env,
        )
        run_cmd(
            [
                sys.executable,
                "-B",
                str(REPO_ROOT / "scripts" / "06_fetch_eqtl_region_hits.py"),
                "--loci",
                str(dirs["processed"] / "x_loci.tsv.gz"),
                "--eqtl-index",
                str(eqtl_index_out),
                "--out",
                str(dirs["interim"] / "eqtl_hits.tsv.gz"),
                "--summary-out",
                str(dirs["interim"] / "eqtl_lookup_summary.tsv"),
            ],
            env=env,
        )
        run_cmd(
            [
                sys.executable,
                "-B",
                str(REPO_ROOT / "scripts" / "07_map_loci_to_genes.py"),
                "--loci",
                str(dirs["processed"] / "x_loci.tsv.gz"),
                "--genes",
                str(genes_out),
                "--eqtl",
                str(dirs["interim"] / "eqtl_hits.tsv.gz"),
                "--eqtl-lookup-summary",
                str(dirs["interim"] / "eqtl_lookup_summary.tsv"),
                "--out-candidates",
                str(dirs["processed"] / "x_gene_candidates.tsv.gz"),
                "--out-summary",
                str(dirs["processed"] / "x_locus_gene_summary.tsv"),
            ],
            env=env,
        )
        run_cmd(
            [
                sys.executable,
                "-B",
                str(REPO_ROOT / "scripts" / "08_export_release_tables.py"),
                "--selected",
                str(selected_path),
                "--loci",
                str(dirs["processed"] / "x_loci.tsv.gz"),
                "--trait-summary",
                str(dirs["processed"] / "x_trait_locus_summary.tsv"),
                "--gene-summary",
                str(dirs["processed"] / "x_locus_gene_summary.tsv"),
                "--gene-candidates",
                str(dirs["processed"] / "x_gene_candidates.tsv.gz"),
                "--outdir",
                str(release_outdir),
            ],
            env=env,
        )
        built_strata.append(stratum)

    if args.skip_frontend_bundle:
        print(f"[done] built strata for {args.panel_id}: {', '.join(built_strata) if built_strata else 'none'}")
        return

    specs = existing_panel_specs(release_root)
    if not specs:
        raise RuntimeError(f"No release directories found under {release_root}; cannot export frontend bundle.")

    export_frontend_bundle(
        specs,
        frontend_outdir,
        min_curated_traits_per_stratum=args.min_curated_traits_per_stratum,
    )
    print(f"[done] built strata for {args.panel_id}: {', '.join(built_strata) if built_strata else 'none'}")
    print(f"[done] frontend bundle: {frontend_outdir}")


if __name__ == "__main__":
    main()
