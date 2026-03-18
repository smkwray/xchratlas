from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd

RSID_RE = re.compile(r"^rs[1-9][0-9]*$", re.IGNORECASE)
NORMALIZED_EQTL_COLUMNS = [
    "study_id",
    "molecular_trait_id",
    "gene_id",
    "variant_id",
    "rsid",
    "pvalue",
]
_NULLISH_BUILD = {"", "na", "n/a", "none", "null", "unknown"}


def normalize_rsid(value: Any) -> str | None:
    """Return canonical lower-case rsID (e.g., rs1234) or None if invalid."""
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    if not RSID_RE.match(text):
        return None
    return f"rs{text[2:]}"


def normalize_rsid_from_locus(
    locus_row: Mapping[str, Any] | pd.Series,
    *,
    candidate_columns: Sequence[str] = ("lead_rsid", "rsid", "lead_variant_rsid"),
) -> str | None:
    for col in candidate_columns:
        if col in locus_row:
            rsid = normalize_rsid(locus_row.get(col))
            if rsid:
                return rsid
    return None


def _normalize_build_name(build: str | None) -> str | None:
    if build is None:
        return None
    text = str(build).strip().lower()
    if text in _NULLISH_BUILD:
        return None
    compact = text.replace("-", "").replace("_", "")
    aliases = {
        "grch37": "GRCh37",
        "hg19": "GRCh37",
        "b37": "GRCh37",
        "grch38": "GRCh38",
        "hg38": "GRCh38",
        "b38": "GRCh38",
    }
    return aliases.get(compact)


def normalize_build_name(build: str | None) -> str | None:
    return _normalize_build_name(build)


def is_region_query_build_safe(
    source_build: str | None,
    target_build: str | None,
    *,
    allow_unknown: bool = False,
) -> bool:
    """A region query is safe only when source and target builds are known and match."""
    source = _normalize_build_name(source_build)
    target = _normalize_build_name(target_build)
    if source is None or target is None:
        return bool(allow_unknown)
    return source == target


def normalized_eqtl_frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    for col in NORMALIZED_EQTL_COLUMNS:
        if col not in frame.columns:
            frame[col] = pd.NA
    return frame[NORMALIZED_EQTL_COLUMNS + [c for c in frame.columns if c not in NORMALIZED_EQTL_COLUMNS]]


def _extract_association_rows(payload: Any) -> list[dict[str, Any]]:
    if payload is None:
        return []
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, Mapping)]
    if not isinstance(payload, Mapping):
        return []

    embedded = payload.get("_embedded")
    if isinstance(embedded, Mapping):
        for key in ("associations", "association"):
            value = embedded.get(key)
            if isinstance(value, list):
                return [row for row in value if isinstance(row, Mapping)]
            if isinstance(value, Mapping):
                return [value]

    for key in ("associations", "association", "results", "data"):
        value = payload.get(key)
        if isinstance(value, list):
            return [row for row in value if isinstance(row, Mapping)]
        if isinstance(value, Mapping):
            return [value]

    if any(k in payload for k in ("study_id", "molecular_trait_id", "gene_id", "variant_id", "pvalue")):
        return [dict(payload)]
    return []


def _first_non_null(*values: Any) -> Any:
    for value in values:
        if value is None:
            continue
        if isinstance(value, float) and math.isnan(value):
            continue
        if isinstance(value, str) and not value.strip():
            continue
        return value
    return None


def parse_eqtl_api_payload(payload: Any, *, lead_rsid: str | None = None) -> pd.DataFrame:
    """
    Parse API JSON payloads into a stable DataFrame schema.
    Handles multiple payload layouts and variant field shapes.
    """
    associations = _extract_association_rows(payload)
    if not associations:
        return normalized_eqtl_frame([])

    rows: list[dict[str, Any]] = []
    normalized_lead = normalize_rsid(lead_rsid)
    for assoc in associations:
        variant = assoc.get("variant")
        if not isinstance(variant, Mapping):
            variant = {}
        gene = assoc.get("gene")
        if not isinstance(gene, Mapping):
            gene = {}
        study = assoc.get("study")
        if not isinstance(study, Mapping):
            study = {}
        molecular_trait = assoc.get("molecular_trait")
        if not isinstance(molecular_trait, Mapping):
            molecular_trait = {}

        rsid = normalize_rsid(
            _first_non_null(
                assoc.get("rsid"),
                assoc.get("snp_id"),
                variant.get("rsid"),
                variant.get("snp_id"),
                assoc.get("snp"),
                normalized_lead,
            )
        )
        row = {
            "study_id": _first_non_null(assoc.get("study_id"), study.get("study_id"), assoc.get("study")),
            "molecular_trait_id": _first_non_null(
                assoc.get("molecular_trait_id"),
                molecular_trait.get("molecular_trait_id"),
                assoc.get("gene_id"),
            ),
            "gene_id": _first_non_null(assoc.get("gene_id"), gene.get("gene_id"), assoc.get("molecular_trait_id")),
            "variant_id": _first_non_null(
                assoc.get("variant_id"),
                variant.get("variant_id"),
                assoc.get("variant"),
                assoc.get("snp"),
            ),
            "rsid": rsid,
            "pvalue": _first_non_null(assoc.get("pvalue"), assoc.get("p_value"), assoc.get("pval")),
        }
        rows.append(row)

    return normalized_eqtl_frame(rows)


def build_locus_query_status(
    *,
    locus_id: str,
    query_mode: str,
    query_status: str,
    query_detail: str = "",
    lead_rsid: str | None = None,
    source_build: str | None = None,
    target_build: str | None = None,
) -> dict[str, Any]:
    return {
        "locus_id": str(locus_id),
        "query_mode": str(query_mode),
        "query_status": str(query_status),
        "query_detail": str(query_detail or ""),
        "lead_rsid": normalize_rsid(lead_rsid),
        "source_build": _normalize_build_name(source_build),
        "target_build": _normalize_build_name(target_build),
    }


def variant_recoder_input_from_locus(locus_row: Mapping[str, Any] | pd.Series) -> str | None:
    chrom = locus_row.get("lead_chr")
    pos = pd.to_numeric(locus_row.get("lead_pos"), errors="coerce")
    ref = locus_row.get("lead_ref")
    alt = locus_row.get("lead_alt")

    if pd.isna(pos) or ref is None or alt is None:
        return None

    chrom_text = str(chrom or "X").strip()
    if not chrom_text:
        chrom_text = "X"
    if chrom_text.lower().startswith("chr"):
        chrom_text = chrom_text[3:]

    ref_text = str(ref).strip()
    alt_text = str(alt).strip()
    if not ref_text or not alt_text:
        return None

    return f"{chrom_text} {int(pos)} . {ref_text} {alt_text}"


def extract_variant_recoder_rsids(payload: Any) -> dict[str, str]:
    resolved: dict[str, str] = {}
    if not isinstance(payload, list):
        return resolved

    for item in payload:
        if not isinstance(item, Mapping):
            continue
        for allele_payload in item.values():
            if not isinstance(allele_payload, Mapping):
                continue
            input_value = allele_payload.get("input")
            if input_value is None:
                continue
            rsids = allele_payload.get("id")
            if not isinstance(rsids, Sequence) or isinstance(rsids, (str, bytes)):
                rsids = [rsids] if rsids else []
            normalized = [normalize_rsid(value) for value in rsids]
            normalized = [value for value in normalized if value]
            if normalized:
                resolved[str(input_value)] = normalized[0]
                break
    return resolved


def eqtl_variant_query_from_components(
    chrom: Any,
    pos: Any,
    ref: Any,
    alt: Any,
) -> str | None:
    pos_num = pd.to_numeric(pos, errors="coerce")
    if pd.isna(pos_num):
        return None

    chrom_text = str(chrom or "").strip()
    ref_text = str(ref or "").strip()
    alt_text = str(alt or "").strip()
    if not chrom_text or not ref_text or not alt_text:
        return None
    if chrom_text.lower().startswith("chr"):
        chrom_text = chrom_text[3:]
    return f"chr{chrom_text}_{int(pos_num)}_{ref_text}_{alt_text}"
