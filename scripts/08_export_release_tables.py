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

from xatlas.io import coerce_bool, write_tsv
from xatlas.scoring import LOOKUP_COMPLETE_STATUSES, compute_trait_score, load_weights, summarize_eqtl_lookup_note
from xatlas.settings import INTERIM_DIR, PROCESSED_DIR, RELEASE_DIR


def _status_count(locus_rows: pd.DataFrame, status: str) -> int:
    if "eqtl_lookup_status" not in locus_rows.columns:
        return 0
    return int(locus_rows["eqtl_lookup_status"].astype(str).str.lower().eq(status).sum())


def build_trait_gene_rollup(gene_candidates: pd.DataFrame) -> pd.DataFrame:
    if gene_candidates.empty:
        return pd.DataFrame()

    ranked = gene_candidates.sort_values(["trait_id", "candidate_rank"]).copy()
    lookup_source = ranked["eqtl_lookup_hit"] if "eqtl_lookup_hit" in ranked.columns else pd.Series([False] * len(ranked), index=ranked.index)
    supported_source = ranked["eqtl_supported"] if "eqtl_supported" in ranked.columns else pd.Series([False] * len(ranked), index=ranked.index)
    ranked["eqtl_lookup_hit_bool"] = lookup_source.map(lambda value: coerce_bool(value) is True)
    ranked["eqtl_supported_bool"] = supported_source.map(lambda value: coerce_bool(value) is True)
    if "eqtl_gene_role" in ranked.columns:
        ranked["eqtl_gene_role"] = ranked["eqtl_gene_role"].fillna("candidate").astype(str)
    else:
        ranked["eqtl_gene_role"] = "candidate"
    top_gene_rows = []
    for trait_id, group in ranked.groupby("trait_id", sort=False):
        candidate_group = group.loc[group["eqtl_gene_role"].str.lower().ne("followup")].copy()
        top_group = candidate_group if not candidate_group.empty else group
        top_row = top_group.iloc[0]
        lookup_note = summarize_eqtl_lookup_note(group)
        locus_rows = group.sort_values(["locus_id", "candidate_rank"]).drop_duplicates(subset=["locus_id"], keep="first")
        completed_mask = locus_rows["eqtl_lookup_status"].astype(str).str.lower().isin(LOOKUP_COMPLETE_STATUSES)
        supported_gene_mask = candidate_group["eqtl_supported_bool"] if not candidate_group.empty else pd.Series(dtype=bool)
        lookup_hit_gene_mask = group["eqtl_lookup_hit_bool"]

        row = {
            "trait_id": trait_id,
            "top_candidate_genes": ",".join(pd.Series(top_group["gene_name"]).dropna().astype(str).head(3)),
            "eqtl_lookup_hit": bool(lookup_hit_gene_mask.any()),
            "eqtl_lookup_hit_gene_count": int(group.loc[lookup_hit_gene_mask, "gene_id"].dropna().astype(str).nunique()) if "gene_id" in group.columns else 0,
            "eqtl_lookup_hit_locus_count": int(group.loc[lookup_hit_gene_mask, "locus_id"].nunique()),
            "eqtl_lookup_hit_count": int(pd.to_numeric(locus_rows["eqtl_lookup_n_hits"], errors="coerce").fillna(0).sum()),
            "eqtl_supported": bool(supported_gene_mask.any()),
            "best_mapping_relation": top_row.get("mapping_relation"),
            "eqtl_supported_gene_count": int(candidate_group.loc[supported_gene_mask, "gene_id"].dropna().astype(str).nunique()) if "gene_id" in candidate_group.columns else 0,
            "eqtl_supported_locus_count": int(
                candidate_group.loc[supported_gene_mask, "locus_id"].nunique()
            ),
            "eqtl_lookup_complete_locus_count": int(completed_mask.sum()),
            "eqtl_lookup_no_hit_locus_count": _status_count(locus_rows, "no_hits"),
            "eqtl_lookup_skipped_no_rsid_locus_count": _status_count(locus_rows, "skipped_no_rsid"),
            "eqtl_lookup_incomplete_locus_count": int(
                (~locus_rows["eqtl_lookup_status"].astype(str).str.lower().isin(LOOKUP_COMPLETE_STATUSES | {"skipped_no_rsid"})).sum()
            ),
            "eqtl_total_hit_count": int(pd.to_numeric(locus_rows["eqtl_lookup_n_hits"], errors="coerce").fillna(0).sum()),
        }
        if lookup_note:
            row["eqtl_lookup_note"] = lookup_note
        top_gene_rows.append(row)
    return pd.DataFrame(top_gene_rows)


def build_release_loci(loci: pd.DataFrame, gene_summary: pd.DataFrame) -> pd.DataFrame:
    if loci.empty:
        return loci
    if gene_summary.empty:
        return loci

    keep_cols = [
        "locus_id",
        "top_gene_id",
        "top_gene_name",
        "best_mapping_relation",
        "eqtl_lookup_hit",
        "eqtl_lookup_hit_gene_count",
        "eqtl_supported",
        "eqtl_supported_gene_count",
        "eqtl_study_count",
        "best_eqtl_pvalue",
        "eqtl_lookup_status",
        "eqtl_lookup_mode",
        "eqtl_lookup_n_hits",
    ]
    usable = gene_summary[[c for c in keep_cols if c in gene_summary.columns]].drop_duplicates(subset=["locus_id"])
    return loci.merge(usable, on="locus_id", how="left")


def main() -> None:
    parser = argparse.ArgumentParser(description="Export release tables for the chrX atlas v1.")
    parser.add_argument("--selected", default=str(INTERIM_DIR / "selected_traits.tsv"))
    parser.add_argument("--loci", default=str(PROCESSED_DIR / "x_loci.tsv.gz"))
    parser.add_argument("--trait-summary", default=str(PROCESSED_DIR / "x_trait_locus_summary.tsv"))
    parser.add_argument("--gene-summary", default=str(PROCESSED_DIR / "x_locus_gene_summary.tsv"))
    parser.add_argument("--gene-candidates", default=str(PROCESSED_DIR / "x_gene_candidates.tsv.gz"))
    parser.add_argument("--weights", default="config/scoring_weights.yml")
    parser.add_argument("--outdir", default=str(RELEASE_DIR))
    args = parser.parse_args()
    outdir = Path(args.outdir)

    selected = pd.read_csv(args.selected, sep="\t")
    loci = pd.read_csv(args.loci, sep="\t", compression="gzip") if Path(args.loci).exists() else pd.DataFrame()
    trait_summary = pd.read_csv(args.trait_summary, sep="\t") if Path(args.trait_summary).exists() else pd.DataFrame()
    gene_summary = pd.read_csv(args.gene_summary, sep="\t") if Path(args.gene_summary).exists() else pd.DataFrame()
    gene_candidates = pd.read_csv(args.gene_candidates, sep="\t", compression="gzip") if Path(args.gene_candidates).exists() else pd.DataFrame()

    weights = load_weights(args.weights)

    top_gene_by_trait = build_trait_gene_rollup(gene_candidates)
    release_loci = build_release_loci(loci, gene_summary)

    score_rows = []
    summary_lookup = trait_summary.set_index("trait_id") if not trait_summary.empty else None
    gene_lookup = top_gene_by_trait.set_index("trait_id") if not top_gene_by_trait.empty else None

    for _, row in selected.iterrows():
        trait_id = row["trait_id"]
        summary_row = summary_lookup.loc[trait_id].to_dict() if summary_lookup is not None and trait_id in summary_lookup.index else {
            "trait_id": trait_id,
            "n_loci": 0,
            "lead_neglog10_pvalue": pd.NA,
            "pvalue_scale": "",
            "par_only_signal": False,
            "any_nonpar_locus": False,
        }
        gene_row = gene_lookup.loc[trait_id].to_dict() if gene_lookup is not None and trait_id in gene_lookup.index else {}
        score = compute_trait_score(row.to_dict(), summary_row, gene_row=gene_row, weights=weights)
        score_rows.append({**row.to_dict(), **summary_row, **gene_row, **score})

    trait_scores = pd.DataFrame(score_rows).sort_values(
        ["x_evidence_score", "coverage_grade", "domain", "description"],
        ascending=[False, True, True, True],
    )
    write_tsv(trait_scores, outdir / "trait_scores.tsv")
    if not release_loci.empty:
        write_tsv(release_loci, outdir / "x_loci.tsv")
    if not gene_candidates.empty:
        write_tsv(gene_candidates, outdir / "x_gene_candidates.tsv")
    write_tsv(selected, outdir / "manifest_snapshot.tsv")

    readme = outdir / "README.txt"
    readme.write_text(
        "\n".join(
            [
                "chrXatlas v1 release tables",
                f"traits: {len(trait_scores):,}",
                f"loci: {len(loci):,}",
                f"gene candidates: {len(gene_candidates):,}",
                "",
                "Primary file: trait_scores.tsv",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"[done] release tables exported to {outdir}/")


if __name__ == "__main__":
    main()
