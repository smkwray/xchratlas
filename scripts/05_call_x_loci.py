#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse
from pathlib import Path

import pandas as pd

from xatlas.io import write_tsv
from xatlas.loci import call_loci_from_dataframe
from xatlas.settings import INTERIM_DIR, PROCESSED_DIR, RAW_DIR
from xatlas.variant_bridge import build_variant_bridge


def main() -> None:
    parser = argparse.ArgumentParser(description="Call simple chrX loci from reduced Pan-UKB extracts.")
    parser.add_argument("--extract-dir", default=str(INTERIM_DIR / "panukb_x"))
    parser.add_argument("--window-bp", type=int, default=500_000)
    parser.add_argument("--selected", default=str(INTERIM_DIR / "selected_traits.tsv"))
    parser.add_argument("--variant-manifest", default=str(RAW_DIR / "panukb" / "full_variant_qc_metrics.txt.bgz"))
    parser.add_argument("--out-loci", default=str(PROCESSED_DIR / "x_loci.tsv.gz"))
    parser.add_argument("--out-trait-summary", default=str(PROCESSED_DIR / "x_trait_locus_summary.tsv"))
    args = parser.parse_args()

    extract_dir = Path(args.extract_dir)
    selected = pd.read_csv(args.selected, sep="\t")
    trait_lookup = selected.set_index("trait_id") if "trait_id" in selected.columns else None

    locus_frames = []
    summaries = []

    trait_ids = list(selected["trait_id"]) if "trait_id" in selected.columns else [
        path.name.replace(".tsv.gz", "") for path in sorted(extract_dir.glob("*.tsv.gz"))
    ]

    for trait_id in trait_ids:
        path = extract_dir / f"{trait_id}.tsv.gz"
        if not path.exists():
            print(f"[warn] {trait_id} skipped: missing extract {path}")
            summaries.append(
                {
                    "trait_id": trait_id,
                    "n_sig_variants": 0,
                    "n_loci": 0,
                    "lead_neglog10_pvalue": pd.NA,
                    "pvalue_column": pd.NA,
                    "pvalue_scale": pd.NA,
                    "par_only_signal": False,
                    "any_nonpar_locus": False,
                    "processing_note": f"missing extract {path}",
                }
            )
            continue
        print(f"[loci] {trait_id}")
        df = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
        try:
            loci_df, summary = call_loci_from_dataframe(df, trait_id=trait_id, window_bp=args.window_bp)
        except Exception as exc:
            print(f"[warn] {trait_id} skipped: {exc}")
            loci_df = pd.DataFrame()
            summary = {
                "trait_id": trait_id,
                "n_sig_variants": 0,
                "n_loci": 0,
                "lead_neglog10_pvalue": pd.NA,
                "pvalue_column": pd.NA,
                "pvalue_scale": pd.NA,
                "par_only_signal": False,
                "any_nonpar_locus": False,
                "processing_note": str(exc),
            }

        if not loci_df.empty and trait_lookup is not None and trait_id in trait_lookup.index:
            meta = trait_lookup.loc[trait_id]
            for col in ["query_id", "domain", "description", "trait_type", "pheno_sex"]:
                if col in meta.index:
                    loci_df[col] = meta[col]

        locus_frames.append(loci_df)
        summaries.append(summary)

    all_loci = pd.concat(locus_frames, ignore_index=True) if locus_frames else pd.DataFrame()
    all_summaries = pd.DataFrame(summaries)

    variant_manifest = Path(args.variant_manifest)
    if not all_loci.empty and variant_manifest.exists():
        print(f"[bridge] annotating lead variants from {variant_manifest}")
        all_loci, bridge_stats = build_variant_bridge(all_loci, variant_manifest)
        print(
            "[bridge] matched "
            f"{bridge_stats['matched_keys']}/{bridge_stats['unique_keys']} unique keys; "
            f"{bridge_stats['matched_rsid_loci']:,}/{len(all_loci):,} loci have true lead_rsid; "
            f"{bridge_stats['matched_loci']:,}/{len(all_loci):,} loci have lead_varid"
        )

    write_tsv(all_loci, args.out_loci)
    write_tsv(all_summaries, args.out_trait_summary)

    print(f"[done] loci: {len(all_loci):,}")
    print(f"[done] trait summaries: {len(all_summaries):,}")


if __name__ == "__main__":
    main()
