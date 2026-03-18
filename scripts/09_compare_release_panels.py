#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse

import pandas as pd

from xatlas.io import coerce_bool, write_tsv
from xatlas.settings import RELEASE_DIR

META_COLUMNS = ["trait_id", "query_id", "domain", "description"]
NUMERIC_COLUMNS = [
    "n_loci",
    "x_evidence_score",
    "eqtl_supported_locus_count",
    "eqtl_total_hit_count",
]
TEXT_COLUMNS = ["coverage_grade"]
BOOL_COLUMNS = ["eqtl_supported"]


def _prepare_release(trait_scores: pd.DataFrame, *, suffix: str, presence_column: str) -> pd.DataFrame:
    cols = [c for c in META_COLUMNS + NUMERIC_COLUMNS + TEXT_COLUMNS + BOOL_COLUMNS if c in trait_scores.columns]
    frame = trait_scores[cols].copy()
    frame[presence_column] = True

    if "eqtl_supported" in frame.columns:
        frame["eqtl_supported"] = frame["eqtl_supported"].map(coerce_bool).fillna(False).astype(bool)

    for col in NUMERIC_COLUMNS:
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")

    rename_map = {col: f"{col}_{suffix}" for col in frame.columns if col not in {"trait_id", presence_column}}
    return frame.rename(columns=rename_map)


def build_panel_release_comparison(small_trait_scores: pd.DataFrame, expanded_trait_scores: pd.DataFrame) -> pd.DataFrame:
    small = _prepare_release(small_trait_scores, suffix="small", presence_column="in_small_release")
    expanded = _prepare_release(expanded_trait_scores, suffix="expanded", presence_column="in_panel_expanded")

    merged = small.merge(expanded, on="trait_id", how="outer")

    merged["query_id"] = merged.get("query_id_small").combine_first(merged.get("query_id_expanded"))
    merged["domain"] = merged.get("domain_small").combine_first(merged.get("domain_expanded"))
    merged["description"] = merged.get("description_small").combine_first(merged.get("description_expanded"))

    for col in ["in_small_release", "in_panel_expanded", "eqtl_supported_small", "eqtl_supported_expanded"]:
        if col in merged.columns:
            merged[col] = merged[col].fillna(False).astype(bool)

    for col in [
        "n_loci_small",
        "n_loci_expanded",
        "x_evidence_score_small",
        "x_evidence_score_expanded",
        "eqtl_supported_locus_count_small",
        "eqtl_supported_locus_count_expanded",
        "eqtl_total_hit_count_small",
        "eqtl_total_hit_count_expanded",
    ]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0)

    for col in ["coverage_grade_small", "coverage_grade_expanded"]:
        if col in merged.columns:
            merged[col] = merged[col].fillna("")

    merged["delta_n_loci"] = merged["n_loci_expanded"] - merged["n_loci_small"]
    merged["delta_x_evidence_score"] = merged["x_evidence_score_expanded"] - merged["x_evidence_score_small"]
    merged["eqtl_supported_delta"] = merged["eqtl_supported_expanded"].astype(int) - merged["eqtl_supported_small"].astype(int)
    merged["delta_eqtl_supported_locus_count"] = (
        merged["eqtl_supported_locus_count_expanded"] - merged["eqtl_supported_locus_count_small"]
    )
    merged["delta_eqtl_total_hit_count"] = merged["eqtl_total_hit_count_expanded"] - merged["eqtl_total_hit_count_small"]

    ordered_cols = [
        "trait_id",
        "query_id",
        "domain",
        "description",
        "in_small_release",
        "in_panel_expanded",
        "n_loci_small",
        "n_loci_expanded",
        "delta_n_loci",
        "x_evidence_score_small",
        "x_evidence_score_expanded",
        "delta_x_evidence_score",
        "eqtl_supported_small",
        "eqtl_supported_expanded",
        "eqtl_supported_delta",
        "eqtl_supported_locus_count_small",
        "eqtl_supported_locus_count_expanded",
        "delta_eqtl_supported_locus_count",
        "eqtl_total_hit_count_small",
        "eqtl_total_hit_count_expanded",
        "delta_eqtl_total_hit_count",
        "coverage_grade_small",
        "coverage_grade_expanded",
    ]

    comparison = merged[[c for c in ordered_cols if c in merged.columns]].copy()
    comparison = comparison.sort_values(
        ["in_panel_expanded", "in_small_release", "domain", "query_id", "trait_id"],
        ascending=[False, False, True, True, True],
        na_position="last",
    ).reset_index(drop=True)
    return comparison


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare the small and expanded chrXatlas release tiers.")
    parser.add_argument("--small-release-dir", default=str(RELEASE_DIR))
    parser.add_argument("--expanded-release-dir", default=str(RELEASE_DIR / "panel_expanded"))
    parser.add_argument("--out", default=str(RELEASE_DIR / "panel_release_comparison.tsv"))
    args = parser.parse_args()

    small_trait_scores = pd.read_csv(Path(args.small_release_dir) / "trait_scores.tsv", sep="\t")
    expanded_trait_scores = pd.read_csv(Path(args.expanded_release_dir) / "trait_scores.tsv", sep="\t")
    comparison = build_panel_release_comparison(small_trait_scores, expanded_trait_scores)
    write_tsv(comparison, args.out)
    print(f"[done] wrote comparison table to {args.out}")


if __name__ == "__main__":
    main()
