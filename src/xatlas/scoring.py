from __future__ import annotations

from typing import Any

import pandas as pd
import yaml

from .io import coerce_bool

DEFAULT_WEIGHTS = {
    "coverage": {"independent_set": 12, "multi_pop_qc": 8, "large_hq_cohort": 10},
    "signal": {"any_locus": 12, "per_extra_locus": 4, "loci_cap": 16, "lead_signal_cap": 20},
    "mapping": {"nearest_gene": 8, "overlapping_gene": 12, "eqtl_support": 15},
    "region": {"nonpar_bonus": 10, "par_only_bonus": 3},
    "penalties": {"no_locus": -10},
    "tiers": {"A": 75, "B": 55, "C": 35},
}

LOOKUP_COMPLETE_STATUSES = {
    "complete",
    "completed",
    "done",
    "found",
    "ok",
    "success",
    "succeeded",
    "not_found",
    "no_hits",
    "empty",
}


def _normalize_token(value: Any) -> str:
    if value is None:
        return ""
    if pd.isna(value):
        return ""
    return str(value).strip().lower()


def eqtl_lookup_note(row: pd.Series | dict[str, Any]) -> str:
    provided_note = row.get("eqtl_lookup_note", "")
    if provided_note is not None and not pd.isna(provided_note):
        provided_note = str(provided_note).strip()
        if provided_note:
            return provided_note

    status = _normalize_token(row.get("eqtl_lookup_status"))
    mode = _normalize_token(row.get("eqtl_lookup_mode"))
    detail = _normalize_token(row.get("eqtl_lookup_detail"))
    n_hits = pd.to_numeric(row.get("eqtl_lookup_n_hits"), errors="coerce")

    has_lookup_data = bool(status or detail or pd.notna(n_hits))
    if not has_lookup_data:
        return ""

    if "provisional" in mode:
        return "eQTL follow-up used provisional region lookup; interpret cautiously."

    if status in LOOKUP_COMPLETE_STATUSES:
        if pd.notna(n_hits) and float(n_hits) == 0.0:
            return "eQTL follow-up: no support observed in lookup summary."
        if pd.notna(n_hits) and float(n_hits) > 0.0 and coerce_bool(row.get("eqtl_supported")) is False:
            return "eQTL follow-up returned lookup-hit associations, but none met the strict candidate-gene support rule."
        if pd.notna(n_hits):
            return ""
        return "eQTL follow-up: not safely assessed / lookup incomplete."

    return "eQTL follow-up: not safely assessed / lookup incomplete."


def summarize_eqtl_lookup_note(rows: pd.DataFrame) -> str:
    if rows.empty:
        return ""
    if not {"eqtl_lookup_status", "eqtl_lookup_detail", "eqtl_lookup_n_hits"}.intersection(rows.columns):
        return ""

    locus_rows = rows.drop_duplicates(subset=["locus_id"]) if "locus_id" in rows.columns else rows.drop_duplicates()

    saw_complete_lookup = False
    all_zero_hits = True
    saw_provisional_mode = False
    saw_skipped_no_rsid = False
    saw_failed_lookup = False
    for _, row in locus_rows.iterrows():
        status = _normalize_token(row.get("eqtl_lookup_status"))
        mode = _normalize_token(row.get("eqtl_lookup_mode"))
        detail = _normalize_token(row.get("eqtl_lookup_detail"))
        n_hits = pd.to_numeric(row.get("eqtl_lookup_n_hits"), errors="coerce")

        if not (status or detail or pd.notna(n_hits)):
            return "eQTL follow-up: not safely assessed / lookup incomplete."
        if "provisional" in mode:
            saw_provisional_mode = True
        if status in LOOKUP_COMPLETE_STATUSES:
            if pd.isna(n_hits):
                return "eQTL follow-up: not safely assessed / lookup incomplete."
            saw_complete_lookup = True
            if float(n_hits) > 0:
                all_zero_hits = False
            continue
        if status == "skipped_no_rsid":
            saw_skipped_no_rsid = True
            continue
        saw_failed_lookup = True

    if saw_provisional_mode:
        return "eQTL follow-up used provisional region lookup; interpret cautiously."
    if saw_complete_lookup and saw_failed_lookup:
        return "eQTL follow-up completed for some loci, but some loci could not be safely assessed."
    if saw_complete_lookup and saw_skipped_no_rsid:
        if all_zero_hits:
            return "eQTL follow-up found no support among rsID-addressable loci; some loci lacked a usable rsID."
        return "eQTL follow-up completed for rsID-addressable loci; some loci lacked a usable rsID."
    if saw_failed_lookup:
        return "eQTL follow-up: not safely assessed / lookup incomplete."
    if saw_complete_lookup and all_zero_hits:
        return "eQTL follow-up: no support observed in lookup summary."
    if saw_complete_lookup and not all_zero_hits and not rows.get("eqtl_supported", pd.Series(dtype=bool)).map(coerce_bool).fillna(False).any():
        return "eQTL follow-up returned lookup-hit associations, but none met the strict candidate-gene support rule."
    return ""


def load_weights(path: str | None = None) -> dict[str, Any]:
    if not path:
        return DEFAULT_WEIGHTS
    with open(path, "rt", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    merged = DEFAULT_WEIGHTS.copy()
    for k, v in (loaded or {}).items():
        if isinstance(v, dict) and isinstance(merged.get(k), dict):
            merged[k] = {**merged[k], **v}
        else:
            merged[k] = v
    return merged


def coverage_grade(row: pd.Series | dict[str, Any]) -> str:
    get = row.get if isinstance(row, dict) else row.__getitem__
    in_independent = coerce_bool(get("in_max_independent_set")) is True
    num_qc = pd.to_numeric(get("num_pops_pass_qc"), errors="coerce")
    hq_n = pd.to_numeric(get("n_cases_hq_cohort_both_sexes"), errors="coerce")

    if in_independent and pd.notna(num_qc) and num_qc >= 2 and pd.notna(hq_n) and hq_n >= 50_000:
        return "A"
    if (in_independent and pd.notna(hq_n) and hq_n >= 10_000) or (pd.notna(num_qc) and num_qc >= 2):
        return "B"
    if pd.notna(hq_n) and hq_n >= 1_000:
        return "C"
    return "U"


def tier_from_score(score: float, weights: dict[str, Any]) -> str:
    tiers = weights["tiers"]
    if score >= tiers["A"]:
        return "A"
    if score >= tiers["B"]:
        return "B"
    if score >= tiers["C"]:
        return "C"
    return "D"


def compute_trait_score(
    manifest_row: pd.Series | dict[str, Any],
    summary_row: pd.Series | dict[str, Any],
    *,
    gene_row: pd.Series | dict[str, Any] | None = None,
    weights: dict[str, Any] | None = None,
) -> dict[str, Any]:
    weights = weights or DEFAULT_WEIGHTS
    m = manifest_row
    s = summary_row
    g = gene_row or {}

    score = 0.0
    notes: list[str] = []

    in_independent = coerce_bool(m.get("in_max_independent_set")) is True
    num_qc = pd.to_numeric(m.get("num_pops_pass_qc"), errors="coerce")
    hq_n = pd.to_numeric(m.get("n_cases_hq_cohort_both_sexes"), errors="coerce")

    if in_independent:
        score += weights["coverage"]["independent_set"]
    if pd.notna(num_qc) and num_qc >= 2:
        score += weights["coverage"]["multi_pop_qc"]
    if pd.notna(hq_n) and hq_n >= 50_000:
        score += weights["coverage"]["large_hq_cohort"]

    n_loci = int(pd.to_numeric(s.get("n_loci"), errors="coerce") or 0)
    lead_score = pd.to_numeric(s.get("lead_neglog10_pvalue"), errors="coerce")

    if n_loci > 0:
        score += weights["signal"]["any_locus"]
        score += min(weights["signal"]["loci_cap"], weights["signal"]["per_extra_locus"] * max(0, n_loci - 1))
        if pd.notna(lead_score):
            score += min(weights["signal"]["lead_signal_cap"], max(0.0, float(lead_score) - 7.30103))
    else:
        score += weights["penalties"]["no_locus"]
        notes.append("No genome-wide significant chrX loci in current build.")

    if bool(s.get("any_nonpar_locus")):
        score += weights["region"]["nonpar_bonus"]
    elif bool(s.get("par_only_signal")) and n_loci > 0:
        score += weights["region"]["par_only_bonus"]
        notes.append("Signal is PAR-only; interpret cautiously.")

    mapping_relation = str(g.get("best_mapping_relation", "") or "")
    eqtl_supported = bool(g.get("eqtl_supported", False))

    if mapping_relation.startswith("lead_overlap") or mapping_relation.startswith("locus_overlap"):
        score += weights["mapping"]["overlapping_gene"]
    elif mapping_relation:
        score += weights["mapping"]["nearest_gene"]

    if eqtl_supported:
        score += weights["mapping"]["eqtl_support"]

    pheno_sex = str(m.get("pheno_sex", "") or "")
    if pheno_sex == "both_sexes":
        notes.append("GWAS is sex-combined; dosage-compensation assumptions are not explicit.")

    pscale = str(s.get("pvalue_scale", "") or "")
    if pscale == "ln":
        notes.append("P-value scale was inferred as ln(P); check whether archived schema was used.")
    elif not pscale:
        notes.append("P-value scale was not recorded.")

    lookup_note = eqtl_lookup_note(g)
    if lookup_note:
        notes.append(lookup_note)

    score = max(0.0, min(100.0, score))
    tier = tier_from_score(score, weights)
    return {
        "coverage_grade": coverage_grade(m),
        "x_evidence_score": round(score, 2),
        "tier": tier,
        "confidence_notes": " ".join(dict.fromkeys(notes)).strip(),
    }
