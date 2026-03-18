#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse

import pandas as pd

from xatlas.io import write_tsv
from xatlas.settings import INTERIM_DIR, PROCESSED_DIR


def strip_version(gene_id: str | None) -> str | None:
    if gene_id is None:
        return None
    return str(gene_id).split(".")[0]


def load_eqtl_lookup_summary(path: str | Path) -> pd.DataFrame | None:
    path = Path(path)
    if not path.exists():
        return None
    try:
        summary = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return None
    if "locus_id" not in summary.columns:
        return None
    summary = summary.copy()
    summary = summary.loc[summary["locus_id"].notna()].copy()
    summary["_locus_id_key"] = summary["locus_id"].map(lambda value: str(value))
    return summary.drop_duplicates("_locus_id_key", keep="last").set_index("_locus_id_key")


def gene_distance_to_lead(gene: pd.Series, lead_pos: int) -> int:
    if int(gene["start"]) <= lead_pos <= int(gene["end"]):
        return 0
    return min(abs(int(gene["start"]) - lead_pos), abs(int(gene["end"]) - lead_pos))


def gene_mapping_relation(
    gene: pd.Series,
    *,
    lead_pos: int,
    locus_start: int,
    locus_end: int,
    flank_bp: int,
) -> str:
    start = int(gene["start"])
    end = int(gene["end"])
    if start <= lead_pos <= end:
        return "lead_overlap"
    if end >= locus_start and start <= locus_end:
        return "locus_overlap"
    if end >= lead_pos - flank_bp and start <= lead_pos + flank_bp:
        return "nearby_window"
    return "eqtl_hit_gene"


def append_missing_eqtl_hit_genes(
    base_subset: pd.DataFrame,
    *,
    genes: pd.DataFrame,
    eqtl_support: dict[str, dict[str, object]],
    lead_pos: int,
    locus_start: int,
    locus_end: int,
    flank_bp: int,
) -> pd.DataFrame:
    if base_subset.empty or not eqtl_support:
        return base_subset

    existing = set(base_subset["gene_id_clean"].dropna().astype(str))
    missing_gene_ids = [gene_id for gene_id in eqtl_support if gene_id and gene_id not in existing]
    if not missing_gene_ids:
        return base_subset

    extras = genes.loc[genes["gene_id_clean"].isin(missing_gene_ids)].copy()
    if extras.empty:
        return base_subset

    extras["distance_to_lead_bp"] = extras.apply(gene_distance_to_lead, axis=1, lead_pos=lead_pos)
    extras["mapping_relation"] = extras.apply(
        gene_mapping_relation,
        axis=1,
        lead_pos=lead_pos,
        locus_start=locus_start,
        locus_end=locus_end,
        flank_bp=flank_bp,
    )
    extras["best_eqtl_pvalue"] = extras["gene_id_clean"].map(
        lambda gene_id: eqtl_support.get(gene_id, {}).get("best_eqtl_pvalue", pd.NA)
    )
    extras = extras.sort_values(["best_eqtl_pvalue", "distance_to_lead_bp", "gene_biotype", "gene_name"], na_position="last")
    extras = extras.drop_duplicates(subset=["gene_id_clean"], keep="first")
    return pd.concat([base_subset, extras], ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Map chrX loci to candidate genes.")
    parser.add_argument("--loci", default=str(PROCESSED_DIR / "x_loci.tsv.gz"))
    parser.add_argument("--genes", default=str(INTERIM_DIR / "ensembl_x_genes.tsv"))
    parser.add_argument("--eqtl", default=str(INTERIM_DIR / "eqtl_hits.tsv.gz"))
    parser.add_argument("--eqtl-lookup-summary", default=str(INTERIM_DIR / "eqtl_lookup_summary.tsv"))
    parser.add_argument("--flank-bp", type=int, default=100_000)
    parser.add_argument("--out-candidates", default=str(PROCESSED_DIR / "x_gene_candidates.tsv.gz"))
    parser.add_argument("--out-summary", default=str(PROCESSED_DIR / "x_locus_gene_summary.tsv"))
    args = parser.parse_args()

    loci = pd.read_csv(args.loci, sep="\t", compression="gzip")
    genes = pd.read_csv(args.genes, sep="\t")

    eqtl = None
    if Path(args.eqtl).exists():
        try:
            eqtl = pd.read_csv(args.eqtl, sep="\t", compression="gzip")
        except pd.errors.EmptyDataError:
            eqtl = None
        if eqtl is not None:
            gene_like = None
            for candidate in ["gene_id", "molecular_trait_id"]:
                if candidate in eqtl.columns:
                    gene_like = candidate
                    break
            if gene_like:
                eqtl["gene_id_clean"] = eqtl[gene_like].map(strip_version)

    eqtl_lookup_summary = load_eqtl_lookup_summary(args.eqtl_lookup_summary)

    rows = []

    genes = genes.copy()
    genes["gene_id_clean"] = genes["gene_id"].map(strip_version)

    for _, locus in loci.iterrows():
        trait_id = locus["trait_id"]
        locus_id = locus["locus_id"]
        lead_pos = int(locus["lead_pos"])
        locus_start = int(locus["locus_start"])
        locus_end = int(locus["locus_end"])

        subset = genes.loc[(genes["end"] >= locus_start) & (genes["start"] <= locus_end)].copy()
        mapping_relation = "locus_overlap"

        if subset.empty:
            subset = genes.loc[
                (genes["end"] >= lead_pos - args.flank_bp) & (genes["start"] <= lead_pos + args.flank_bp)
            ].copy()
            mapping_relation = "nearby_window"

        if subset.empty:
            fallback = genes.copy()
            fallback["distance_to_lead_bp"] = fallback.apply(
                lambda r: min(abs(int(r["start"]) - lead_pos), abs(int(r["end"]) - lead_pos)),
                axis=1,
            )
            subset = fallback.nsmallest(5, "distance_to_lead_bp").copy()
            mapping_relation = "nearest_gene"

        subset["distance_to_lead_bp"] = subset.apply(
            lambda r: 0 if int(r["start"]) <= lead_pos <= int(r["end"]) else min(abs(int(r["start"]) - lead_pos), abs(int(r["end"]) - lead_pos)),
            axis=1,
        )
        subset["mapping_relation"] = subset.apply(
            lambda r: "lead_overlap" if int(r["start"]) <= lead_pos <= int(r["end"]) else mapping_relation,
            axis=1,
        )

        eqtl_support = {}
        if eqtl is not None and not eqtl.empty:
            locus_eqtl = eqtl.loc[eqtl["locus_id"] == locus_id].copy()
            if "gene_id_clean" in locus_eqtl.columns:
                for gene_id, block in locus_eqtl.groupby("gene_id_clean"):
                    pcol = "pvalue" if "pvalue" in block.columns else None
                    best_p = pd.to_numeric(block[pcol], errors="coerce").min() if pcol else pd.NA
                    eqtl_support[gene_id] = {
                        "eqtl_supported": True,
                        "eqtl_study_count": int(block["study_id"].nunique()) if "study_id" in block.columns else pd.NA,
                        "best_eqtl_pvalue": best_p,
                    }

        lookup_info = {
            "eqtl_lookup_status": pd.NA,
            "eqtl_lookup_mode": pd.NA,
            "eqtl_lookup_detail": pd.NA,
            "eqtl_lookup_n_hits": pd.NA,
        }
        if eqtl_lookup_summary is not None:
            key = str(locus_id)
            if key in eqtl_lookup_summary.index:
                lookup_row = eqtl_lookup_summary.loc[key]
                if isinstance(lookup_row, pd.DataFrame):
                    lookup_row = lookup_row.iloc[-1]
                lookup_info = {
                    "eqtl_lookup_status": lookup_row.get("query_status", pd.NA),
                    "eqtl_lookup_mode": lookup_row.get("query_mode", pd.NA),
                    "eqtl_lookup_detail": lookup_row.get("query_detail", pd.NA),
                    "eqtl_lookup_n_hits": pd.to_numeric(lookup_row.get("n_hits"), errors="coerce"),
                }

        subset = subset.sort_values(["distance_to_lead_bp", "gene_biotype", "gene_name"])
        base_subset = subset.head(5).copy()
        candidate_subset = append_missing_eqtl_hit_genes(
            base_subset,
            genes=genes,
            eqtl_support=eqtl_support,
            lead_pos=lead_pos,
            locus_start=locus_start,
            locus_end=locus_end,
            flank_bp=args.flank_bp,
        )

        for rank, (_, gene) in enumerate(candidate_subset.iterrows(), start=1):
            info = eqtl_support.get(gene["gene_id_clean"], {})
            rows.append(
                {
                    "trait_id": trait_id,
                    "locus_id": locus_id,
                    "candidate_rank": rank,
                    "gene_id": gene["gene_id"],
                    "gene_name": gene.get("gene_name"),
                    "gene_biotype": gene.get("gene_biotype"),
                    "mapping_relation": gene["mapping_relation"],
                    "distance_to_lead_bp": int(gene["distance_to_lead_bp"]),
                    "eqtl_supported": bool(info.get("eqtl_supported", False)),
                    "eqtl_study_count": info.get("eqtl_study_count", 0),
                    "best_eqtl_pvalue": info.get("best_eqtl_pvalue", pd.NA),
                    **lookup_info,
                }
            )

    gene_candidates = pd.DataFrame(rows)
    write_tsv(gene_candidates, args.out_candidates)

    top_by_locus = gene_candidates.sort_values(["locus_id", "candidate_rank"]).groupby("locus_id", as_index=False).first()
    summary = top_by_locus.rename(
        columns={
            "mapping_relation": "best_mapping_relation",
            "gene_name": "top_gene_name",
            "gene_id": "top_gene_id",
        }
    )
    write_tsv(summary, args.out_summary)
    print(f"[done] wrote {len(gene_candidates):,} locus-gene rows")


if __name__ == "__main__":
    main()
