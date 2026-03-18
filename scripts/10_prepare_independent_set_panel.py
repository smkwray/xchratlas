#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse

import numpy as np
import pandas as pd

from xatlas.io import coerce_bool, ensure_dir, maybe_read_bgz_tsv, write_tsv
from xatlas.settings import INTERIM_DIR, RAW_DIR


def make_trait_id(row: pd.Series) -> str:
    pieces = [
        str(row.get("trait_type", "") or ""),
        str(row.get("phenocode", "") or ""),
        str(row.get("pheno_sex", "") or ""),
        str(row.get("coding", "") or ""),
        str(row.get("modifier", "") or ""),
    ]
    cleaned = [re.sub(r"[^A-Za-z0-9._-]+", "-", p.strip()) if p else "na" for p in pieces]
    return "__".join(cleaned)


def normalize_label(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def derive_domain(row: pd.Series) -> str:
    category = str(row.get("category", "") or "")
    if category and ">" in category:
        leaf = category.split(">")[-1].strip()
        if leaf:
            return normalize_label(leaf)
    trait_type = str(row.get("trait_type", "") or "")
    if trait_type:
        return normalize_label(trait_type)
    return "independent_set"


def derive_query_id(row: pd.Series) -> str:
    desc = str(row.get("description", "") or "")
    if desc:
        return normalize_label(desc)[:64] or make_trait_id(row)
    return make_trait_id(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare a selected-traits table for the Pan-UKB max independent set.")
    parser.add_argument("--manifest", default=str(RAW_DIR / "panukb" / "phenotype_manifest.tsv.bgz"))
    parser.add_argument("--outdir", default=str(INTERIM_DIR / "panel_independent_set"))
    parser.add_argument("--min-pops-pass-qc", type=int, default=2)
    parser.add_argument("--require-both-sexes", action="store_true")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    manifest = maybe_read_bgz_tsv(args.manifest)

    independent = manifest["in_max_independent_set"].map(coerce_bool).fillna(False).astype(bool)
    num_pops = pd.to_numeric(manifest.get("num_pops_pass_qc"), errors="coerce").fillna(0)
    mask = independent & (num_pops >= args.min_pops_pass_qc)
    if args.require_both_sexes:
        mask &= manifest.get("pheno_sex", "").fillna("").eq("both_sexes")

    selected = manifest.loc[mask].copy()
    selected["trait_id"] = selected.apply(make_trait_id, axis=1)
    selected["query_id"] = selected.apply(derive_query_id, axis=1)
    selected["domain"] = selected.apply(derive_domain, axis=1)
    selected["priority"] = 1
    selected["query"] = selected["description"].fillna(selected["trait_id"])
    selected["seed_notes"] = "Pan-UKB max independent set"

    hq_cases = pd.to_numeric(selected.get("n_cases_hq_cohort_both_sexes"), errors="coerce").fillna(0)
    selected["selection_score"] = num_pops.loc[selected.index].astype(float) * 100 + np.log10(hq_cases + 1) * 10
    selected = selected.sort_values(
        ["selection_score", "trait_type", "description"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    selected["selection_rank"] = range(1, len(selected) + 1)

    write_tsv(selected, outdir / "selected_traits.tsv")
    print(f"[done] wrote {len(selected):,} independent-set traits to {outdir / 'selected_traits.tsv'}")


if __name__ == "__main__":
    main()
