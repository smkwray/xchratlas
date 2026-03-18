#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse

import pandas as pd

from xatlas.io import ensure_dir, write_tsv


def maybe_concat(paths: list[Path], *, compression: str | None = None) -> pd.DataFrame:
    frames = []
    for path in paths:
        if path.exists():
            frames.append(pd.read_csv(path, sep="\t", compression=compression))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge shard outputs from a parallel chrXatlas panel build.")
    parser.add_argument("--shard-manifest", required=True)
    parser.add_argument("--panel-root", required=True)
    parser.add_argument("--release-outdir", required=True)
    args = parser.parse_args()

    shard_manifest = pd.read_csv(args.shard_manifest, sep="\t")
    panel_root = ensure_dir(args.panel_root)
    release_outdir = ensure_dir(args.release_outdir)

    selected_paths = [Path(p) for p in shard_manifest["selected_path"]]
    shard_dirs = [p.parent for p in selected_paths]

    merged_selected = maybe_concat(selected_paths)
    merged_loci = maybe_concat([d / "processed" / "x_loci.tsv.gz" for d in shard_dirs], compression="gzip")
    merged_trait_summary = maybe_concat([d / "processed" / "x_trait_locus_summary.tsv" for d in shard_dirs])
    merged_eqtl_hits = maybe_concat([d / "interim" / "eqtl_hits.tsv.gz" for d in shard_dirs], compression="gzip")
    merged_eqtl_summary = maybe_concat([d / "interim" / "eqtl_lookup_summary.tsv" for d in shard_dirs])
    merged_gene_candidates = maybe_concat([d / "processed" / "x_gene_candidates.tsv.gz" for d in shard_dirs], compression="gzip")
    merged_gene_summary = maybe_concat([d / "processed" / "x_locus_gene_summary.tsv" for d in shard_dirs])

    interim_dir = ensure_dir(panel_root / "interim")
    processed_dir = ensure_dir(panel_root / "processed")

    write_tsv(merged_selected, interim_dir / "selected_traits.tsv")
    write_tsv(merged_eqtl_hits, interim_dir / "eqtl_hits.tsv.gz")
    write_tsv(merged_eqtl_summary, interim_dir / "eqtl_lookup_summary.tsv")
    write_tsv(merged_loci, processed_dir / "x_loci.tsv.gz")
    write_tsv(merged_trait_summary, processed_dir / "x_trait_locus_summary.tsv")
    write_tsv(merged_gene_candidates, processed_dir / "x_gene_candidates.tsv.gz")
    write_tsv(merged_gene_summary, processed_dir / "x_locus_gene_summary.tsv")

    export_cmd = [
        sys.executable,
        "-B",
        str(REPO_ROOT / "scripts" / "08_export_release_tables.py"),
        "--selected",
        str(interim_dir / "selected_traits.tsv"),
        "--loci",
        str(processed_dir / "x_loci.tsv.gz"),
        "--trait-summary",
        str(processed_dir / "x_trait_locus_summary.tsv"),
        "--gene-summary",
        str(processed_dir / "x_locus_gene_summary.tsv"),
        "--gene-candidates",
        str(processed_dir / "x_gene_candidates.tsv.gz"),
        "--outdir",
        str(release_outdir),
    ]
    subprocess.run(export_cmd, check=True)
    print(f"[done] merged shard outputs into {panel_root} and {release_outdir}")


if __name__ == "__main__":
    main()
