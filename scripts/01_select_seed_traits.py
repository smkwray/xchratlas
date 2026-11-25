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
from xatlas.scoring import effective_hq_case_count
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


def requested_stratum_rows(
    scored: pd.DataFrame,
    requested_pheno_sex: str | None,
    *,
    strict: bool,
    allow_fallback_both_sexes: bool,
) -> tuple[pd.DataFrame, bool]:
    if not requested_pheno_sex:
        return scored, False

    requested = str(requested_pheno_sex or "").strip()
    if not requested:
        return scored, False

    exact = scored.loc[scored["pheno_sex"].fillna("").eq(requested)].copy()
    if requested == "both_sexes" or strict:
        return exact, False
    if not allow_fallback_both_sexes:
        return exact, False

    if not exact.empty:
        return exact, False
    fallback = scored.loc[scored["pheno_sex"].fillna("").eq("both_sexes")].copy()
    return fallback, not fallback.empty


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

    hq_n = effective_hq_case_count(row)
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
    parser.add_argument("--pheno-sex", choices=["both_sexes", "females", "males"], default=None)
    parser.add_argument("--strict-pheno-sex", action="store_true")
    parser.add_argument("--allow-fallback-both-sexes", action="store_true")
    parser.add_argument("--selection-audit-out", default=None)
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    manifest = maybe_read_bgz_tsv(args.manifest)
    panel = pd.read_csv(args.panel)

    manifest["_search_blob"] = build_search_blob(manifest, SEARCH_COLUMNS)

    candidate_rows = []
    selected_rows = []
    audit_rows = []

    for _, seed in panel.sort_values(["priority", "query_id"]).iterrows():
        qid = seed["query_id"]
        query = seed["query"]
        preferred_trait_type = seed.get("preferred_trait_type")
        preferred_pheno_sex = args.pheno_sex or seed.get("preferred_pheno_sex")

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
            audit_rows.append(
                {
                    "query_id": qid,
                    "requested_pheno_sex": preferred_pheno_sex,
                    "selected_trait_id": None,
                    "selected_pheno_sex": None,
                    "selected_description": None,
                    "selection_score": None,
                    "used_fallback": False,
                    "n_candidate_rows_exact_stratum": 0,
                    "n_candidate_rows_any_stratum": 0,
                }
            )
            continue

        n_any_stratum = int(len(scored))
        requested_pool, used_fallback = requested_stratum_rows(
            scored,
            preferred_pheno_sex,
            strict=args.strict_pheno_sex,
            allow_fallback_both_sexes=args.allow_fallback_both_sexes,
        )
        n_exact = 0
        if preferred_pheno_sex:
            n_exact = int(scored["pheno_sex"].fillna("").eq(preferred_pheno_sex).sum())
            scored = requested_pool
            if scored.empty:
                audit_rows.append(
                    {
                        "query_id": qid,
                        "requested_pheno_sex": preferred_pheno_sex,
                        "selected_trait_id": None,
                        "selected_pheno_sex": None,
                        "selected_description": None,
                        "selection_score": None,
                        "used_fallback": used_fallback,
                        "n_candidate_rows_exact_stratum": n_exact,
                        "n_candidate_rows_any_stratum": n_any_stratum,
                    }
                )
                continue

        scored["trait_id"] = scored.apply(make_trait_id, axis=1)
        scored["requested_pheno_sex"] = preferred_pheno_sex
        scored["selection_used_fallback"] = used_fallback
        scored["effective_hq_case_count"] = scored.apply(effective_hq_case_count, axis=1)
        scored = scored.sort_values(
            ["selection_score", "in_max_independent_set", "num_pops_pass_qc", "effective_hq_case_count"],
            ascending=[False, False, False, False],
        )

        topk = scored.head(args.max_per_query).copy()
        topk["selection_rank"] = range(1, len(topk) + 1)
        candidate_rows.append(topk)

        chosen = topk.iloc[0].copy()
        chosen["selection_rank"] = 1
        selected_rows.append(chosen)
        audit_rows.append(
            {
                "query_id": qid,
                "requested_pheno_sex": preferred_pheno_sex,
                "selected_trait_id": chosen.get("trait_id"),
                "selected_pheno_sex": chosen.get("pheno_sex"),
                "selected_description": chosen.get("description"),
                "selection_score": chosen.get("selection_score"),
                "used_fallback": bool(used_fallback),
                "n_candidate_rows_exact_stratum": n_exact,
                "n_candidate_rows_any_stratum": n_any_stratum,
            }
        )

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
    if args.selection_audit_out:
        write_tsv(pd.DataFrame(audit_rows), args.selection_audit_out)

    print(f"[done] candidate rows: {len(candidates):,}")
    print(f"[done] selected rows: {len(selected):,}")


if __name__ == "__main__":
    main()
