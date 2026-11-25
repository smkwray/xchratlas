#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

STRATA = ("both_sexes", "females", "males")
SEARCH_COLUMNS = [
    "description",
    "description_more",
    "category",
    "coding_description",
    "phenocode",
    "trait_type",
    "modifier",
]


def normalize_text(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip().lower())


def tokenize(value: object) -> list[str]:
    return [token for token in re.split(r"[^a-z0-9]+", normalize_text(value)) if token]


def maybe_read_bgz_tsv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    compression = "gzip" if path.suffix.lower() in {".gz", ".bgz"} else "infer"
    return pd.read_csv(path, sep="\t", compression=compression)


def build_search_blob(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    usable = [column for column in columns if column in df.columns]
    if not usable:
        return pd.Series([""] * len(df), index=df.index)
    return df[usable].fillna("").astype(str).agg(" | ".join, axis=1).str.lower()


def effective_hq_cases(row: pd.Series) -> float:
    sex = normalize_text(row.get("pheno_sex"))
    if sex == "females":
        return pd.to_numeric(row.get("n_cases_hq_cohort_females"), errors="coerce")
    if sex == "males":
        return pd.to_numeric(row.get("n_cases_hq_cohort_males"), errors="coerce")
    return pd.to_numeric(row.get("n_cases_hq_cohort_both_sexes"), errors="coerce")


def query_score(row: pd.Series, query: str, preferred_trait_type: str | None, requested_pheno_sex: str) -> float:
    blob = row["_search_blob"]
    desc_norm = normalize_text(row.get("description"))
    query_norm = normalize_text(query)
    query_tokens = tokenize(query)
    score = 0.0

    if desc_norm == query_norm:
        score += 80
    elif desc_norm.startswith(query_norm):
        score += 65
    elif query_norm in desc_norm:
        score += 50
    if query_norm in blob:
        score += 30
    score += 8 * sum(token in blob for token in query_tokens)

    if str(row.get("in_max_independent_set", "")).strip().lower() in {"true", "1", "yes", "y", "t"}:
        score += 20

    num_qc = pd.to_numeric(row.get("num_pops_pass_qc"), errors="coerce")
    if pd.notna(num_qc):
        score += min(float(num_qc), 6.0)

    hq_n = effective_hq_cases(row)
    if pd.notna(hq_n) and float(hq_n) > 0:
        score += min(15.0, np.log10(float(hq_n) + 1) * 3)

    trait_type = str(row.get("trait_type", "") or "")
    if preferred_trait_type and trait_type == preferred_trait_type:
        score += 10

    pheno_sex = str(row.get("pheno_sex", "") or "")
    if pheno_sex == requested_pheno_sex:
        score += 12
    else:
        score -= 100

    modifier = str(row.get("modifier", "") or "")
    if modifier and modifier.lower() not in {"irnt", "na"}:
        score -= 2

    return score


def match_flags(row: pd.Series, query: str, preferred_trait_type: str | None) -> dict[str, object]:
    description = normalize_text(row.get("description"))
    query_norm = normalize_text(query)
    description_tokens = set(tokenize(description))
    query_tokens = set(tokenize(query_norm))

    exact = description == query_norm
    contains = bool(query_norm and (query_norm in description or description in query_norm))
    token_overlap = len(description_tokens & query_tokens)
    token_recall = token_overlap / max(1, len(query_tokens))
    trait_type_match = not preferred_trait_type or str(row.get("trait_type", "") or "") == str(preferred_trait_type)
    curated_match = trait_type_match and (exact or contains or token_recall >= 0.67)

    return {
        "description_exact_match": exact,
        "description_contains_query": contains,
        "query_token_overlap": token_overlap,
        "query_token_recall": round(token_recall, 3),
        "trait_type_match": trait_type_match,
        "curated_match": curated_match,
    }


def collect_stratum_candidates(
    manifest: pd.DataFrame,
    *,
    query: str,
    preferred_trait_type: str | None,
    requested_pheno_sex: str,
    max_per_query: int,
) -> pd.DataFrame:
    scored = manifest.copy()
    scored["selection_score"] = scored.apply(
        query_score,
        axis=1,
        args=(query, preferred_trait_type, requested_pheno_sex),
    )
    scored = scored.loc[scored["selection_score"] >= 1].copy()
    scored = scored.loc[scored["pheno_sex"].fillna("").eq(requested_pheno_sex)].copy()
    if scored.empty:
        return scored

    flags = scored.apply(lambda row: pd.Series(match_flags(row, query, preferred_trait_type)), axis=1)
    scored = pd.concat([scored.reset_index(drop=True), flags.reset_index(drop=True)], axis=1)
    scored["effective_hq_case_count"] = scored.apply(effective_hq_cases, axis=1)
    scored = scored.sort_values(
        [
            "curated_match",
            "description_exact_match",
            "trait_type_match",
            "query_token_recall",
            "selection_score",
            "num_pops_pass_qc",
            "effective_hq_case_count",
        ],
        ascending=[False, False, False, False, False, False, False],
        na_position="last",
    )
    return scored.head(max_per_query).copy()


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit strict sex-stratum coverage for a curated panel config.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--max-per-query", type=int, default=5)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = maybe_read_bgz_tsv(args.manifest)
    panel = pd.read_csv(args.panel)
    manifest["_search_blob"] = build_search_blob(manifest, SEARCH_COLUMNS)

    summary_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []

    for _, seed in panel.sort_values(["priority", "query_id"]).iterrows():
        summary = {
            "query_id": seed["query_id"],
            "domain": seed.get("domain"),
            "query": seed["query"],
            "preferred_trait_type": seed.get("preferred_trait_type"),
        }
        for stratum in STRATA:
            topk = collect_stratum_candidates(
                manifest,
                query=seed["query"],
                preferred_trait_type=seed.get("preferred_trait_type"),
                requested_pheno_sex=stratum,
                max_per_query=args.max_per_query,
            )
            curated = topk.loc[topk["curated_match"]].copy() if not topk.empty else topk
            summary[f"{stratum}_candidate_count"] = int(len(topk))
            summary[f"{stratum}_curated_candidate_count"] = int(len(curated))
            summary[f"{stratum}_available"] = not curated.empty
            if not topk.empty:
                best = topk.iloc[0]
                summary[f"{stratum}_top_description"] = best.get("description")
                summary[f"{stratum}_top_trait_type"] = best.get("trait_type")
                summary[f"{stratum}_top_curated_match"] = bool(best.get("curated_match"))
                summary[f"{stratum}_top_selection_score"] = float(best.get("selection_score"))
                summary[f"{stratum}_top_token_recall"] = float(best.get("query_token_recall"))
                rows = topk.copy()
                rows["requested_pheno_sex"] = stratum
                rows["query_id"] = seed["query_id"]
                candidate_rows.extend(rows.to_dict("records"))
            else:
                summary[f"{stratum}_top_description"] = None
                summary[f"{stratum}_top_trait_type"] = None
                summary[f"{stratum}_top_curated_match"] = False
                summary[f"{stratum}_top_selection_score"] = None
                summary[f"{stratum}_top_token_recall"] = None
        summary_rows.append(summary)

    summary_df = pd.DataFrame(summary_rows)
    candidates_df = pd.DataFrame(candidate_rows)
    summary_df.to_csv(outdir / "stratum_availability_summary.tsv", sep="\t", index=False)
    candidates_df.to_csv(outdir / "stratum_availability_candidates.tsv", sep="\t", index=False)

    print(f"[done] wrote {outdir / 'stratum_availability_summary.tsv'}")
    print(f"[done] wrote {outdir / 'stratum_availability_candidates.tsv'}")
    print()
    print("Curated availability by query:")
    for _, row in summary_df.iterrows():
        statuses = [f"{stratum}={'yes' if bool(row[f'{stratum}_available']) else 'no'}" for stratum in STRATA]
        print(f"  - {row['query_id']}: {', '.join(statuses)}")


if __name__ == "__main__":
    main()
