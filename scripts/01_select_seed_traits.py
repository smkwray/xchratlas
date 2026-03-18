#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from xatlas.io import build_search_blob, coerce_bool, ensure_dir, maybe_read_bgz_tsv, write_tsv
from xatlas.settings import INTERIM_DIR, RAW_DIR

SEARCH_COLUMNS = [
    "description",
    "description_more",
    "category",
    "coding_description",
    "phenocode",
    "trait_type",
    "modifier",
]


def normalize_text(s: str) -> str:
    return re.sub(r"\s+", " ", s.strip().lower())


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


def query_score(row: pd.Series, query: str, preferred_trait_type: str | None, preferred_pheno_sex: str | None) -> float:
    blob = row["_search_blob"]
    desc = str(row.get("description", "") or "").lower()

    q = normalize_text(query)
    desc_norm = normalize_text(str(row.get("description", "") or ""))
    q_tokens = [t for t in re.split(r"[^a-z0-9]+", q) if t]
    token_hits = sum(token in blob for token in q_tokens)

    score = 0.0
    if desc_norm == q:
        score += 80
    elif desc_norm.startswith(q):
        score += 65
    elif q in desc:
        score += 50
    if q in blob:
        score += 30
    score += 8 * token_hits

    if coerce_bool(row.get("in_max_independent_set")) is True:
        score += 20

    num_qc = pd.to_numeric(row.get("num_pops_pass_qc"), errors="coerce")
    if pd.notna(num_qc):
        score += min(float(num_qc), 6.0)

    hq_n = pd.to_numeric(row.get("n_cases_hq_cohort_both_sexes"), errors="coerce")
    if pd.notna(hq_n) and hq_n > 0:
        score += min(15.0, np.log10(float(hq_n) + 1) * 3)

    trait_type = str(row.get("trait_type", "") or "")
    if preferred_trait_type and trait_type == preferred_trait_type:
        score += 10

    pheno_sex = str(row.get("pheno_sex", "") or "")
    if preferred_pheno_sex:
        if pheno_sex == preferred_pheno_sex:
            score += 10
        elif pheno_sex == "both_sexes" and preferred_pheno_sex == "both_sexes":
            score += 8

    modifier = str(row.get("modifier", "") or "")
    if modifier and modifier.lower() not in {"irnt", "na"}:
        score -= 2

    return score


def main() -> None:
    parser = argparse.ArgumentParser(description="Select a balanced small trait panel from the Pan-UKB manifest.")
    parser.add_argument("--manifest", default=str(RAW_DIR / "panukb" / "phenotype_manifest.tsv.bgz"))
    parser.add_argument("--panel", default="config/panel_small.csv")
    parser.add_argument("--outdir", default=str(INTERIM_DIR))
    parser.add_argument("--max-per-query", type=int, default=5)
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    manifest = maybe_read_bgz_tsv(args.manifest)
    panel = pd.read_csv(args.panel)

    manifest["_search_blob"] = build_search_blob(manifest, SEARCH_COLUMNS)

    candidate_rows = []
    selected_rows = []

    for _, seed in panel.sort_values(["priority", "query_id"]).iterrows():
        qid = seed["query_id"]
        query = seed["query"]
        preferred_trait_type = seed.get("preferred_trait_type")
        preferred_pheno_sex = seed.get("preferred_pheno_sex")

        scored = manifest.copy()
        scored["query_id"] = qid
        scored["domain"] = seed["domain"]
        scored["priority"] = seed["priority"]
        scored["query"] = query
        scored["seed_notes"] = seed.get("notes", "")
        scored["selection_score"] = scored.apply(
            query_score,
            axis=1,
            args=(query, preferred_trait_type, preferred_pheno_sex),
        )

        scored = scored.loc[scored["selection_score"] > 0].copy()
        if scored.empty:
            continue

        scored["trait_id"] = scored.apply(make_trait_id, axis=1)
        scored = scored.sort_values(
            ["selection_score", "in_max_independent_set", "num_pops_pass_qc", "n_cases_hq_cohort_both_sexes"],
            ascending=[False, False, False, False],
        )

        topk = scored.head(args.max_per_query).copy()
        topk["selection_rank"] = range(1, len(topk) + 1)
        candidate_rows.append(topk)

        chosen = topk.iloc[0].copy()
        chosen["selection_rank"] = 1
        selected_rows.append(chosen)

    if candidate_rows:
        candidates = pd.concat(candidate_rows, ignore_index=True)
    else:
        candidates = pd.DataFrame()

    if selected_rows:
        selected = pd.DataFrame(selected_rows)
    else:
        selected = pd.DataFrame()

    if not selected.empty:
        # De-duplicate identical trait_ids picked by multiple seed queries, keeping the higher score.
        selected = selected.sort_values(["selection_score", "priority"], ascending=[False, True])
        selected = selected.drop_duplicates(subset=["trait_id"], keep="first").copy()
        selected = selected.sort_values(["domain", "query_id"]).reset_index(drop=True)

    write_tsv(candidates, outdir / "candidate_traits.tsv")
    write_tsv(selected, outdir / "selected_traits.tsv")

    print(f"[done] candidate rows: {len(candidates):,}")
    print(f"[done] selected rows: {len(selected):,}")


if __name__ == "__main__":
    main()
