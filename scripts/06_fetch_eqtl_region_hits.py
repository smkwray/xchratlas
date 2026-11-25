#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse
from typing import Any

import requests
import subprocess
import time
from pathlib import Path

import pandas as pd

from xatlas.constants import (
    ENSEMBL_GRCH37_REST_BASE_URL,
    ENSEMBL_GRCH37_VARIANT_RECODER_URL,
    ENSEMBL_GRCH38_SEQUENCE_BASE_URL,
    EQTL_API_BASE_URL,
)
from xatlas.eqtl import (
    build_locus_query_status,
    eqtl_variant_query_from_components,
    extract_variant_recoder_rsids,
    is_region_query_build_safe,
    normalize_build_name,
    normalize_rsid_from_locus,
    parse_eqtl_api_payload,
    variant_recoder_input_from_locus,
)
from xatlas.io import choose_column, write_tsv
from xatlas.settings import INTERIM_DIR, PROCESSED_DIR


def run_tabix_region(path_or_url: str, region: str) -> pd.DataFrame:
    cmd = ["tabix", "-h", path_or_url, region]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(result.stderr[:1000])
    if not result.stdout.strip():
        return pd.DataFrame()

    lines = result.stdout.splitlines()
    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:] if line.strip()]
    if not rows:
        return pd.DataFrame(columns=header)
    return pd.DataFrame(rows, columns=header)


def compact_detail(text: str, *, limit: int = 240) -> str:
    cleaned = " ".join(str(text or "").split())
    if len(cleaned) <= limit:
        return cleaned
    return cleaned[: limit - 3] + "..."


def extract_message(payload: Any) -> str:
    if isinstance(payload, dict):
        for key in ("message", "detail", "error", "title"):
            value = payload.get(key)
            if value:
                return str(value)
    return ""


def payload_has_association_shape(payload: Any) -> bool:
    if isinstance(payload, list):
        return True
    if not isinstance(payload, dict):
        return False
    embedded = payload.get("_embedded")
    if isinstance(embedded, dict) and any(key in embedded for key in ("associations", "association")):
        return True
    return any(key in payload for key in ("associations", "association", "results", "data"))


def choose_eqtl_study_specs(eqtl_index: pd.DataFrame, max_studies_per_locus: int) -> list[dict[str, str]]:
    if eqtl_index.empty or max_studies_per_locus <= 0:
        return []

    study_col = choose_column(list(eqtl_index.columns), exact=["study_id", "study"], contains_any=["study"])
    dataset_col = choose_column(list(eqtl_index.columns), exact=["dataset_id"], contains_all=["dataset", "id"])
    group_col = choose_column(list(eqtl_index.columns), exact=["qtl_group", "sample_group"], contains_any=["qtl_group", "sample_group"])
    quant_col = choose_column(list(eqtl_index.columns), exact=["quant_method"], contains_all=["quant", "method"])

    if dataset_col is None:
        return []

    eligible = eqtl_index.copy()
    sort_cols = [col for col in ("priority", study_col, dataset_col) if col and col in eligible.columns]
    if sort_cols:
        eligible = eligible.sort_values(sort_cols, na_position="last")

    specs: list[dict[str, str]] = []
    seen: set[str] = set()
    seen_studies: set[str] = set()
    eligible = eligible.loc[eligible[dataset_col].notna()].copy()
    if quant_col:
        ge_mask = eligible[quant_col].astype(str).str.strip().str.lower().eq("ge")
        if ge_mask.any():
            eligible = eligible.loc[ge_mask].copy()

    def maybe_append(row: pd.Series) -> None:
        dataset_id = str(row.get(dataset_col, "") or "").strip()
        if not dataset_id or dataset_id in seen or len(specs) >= max_studies_per_locus:
            return
        study = str(row.get(study_col, "") or "").strip() if study_col else ""
        group = str(row.get(group_col, "") or "").strip() if group_col else ""
        quant = str(row.get(quant_col, "") or "").strip() if quant_col else ""
        specs.append({"dataset_id": dataset_id, "study": study, "qtl_group": group, "quant_method": quant})
        seen.add(dataset_id)
        if study:
            seen_studies.add(study)

    for _, row in eligible.iterrows():
        study = str(row.get(study_col, "") or "").strip() if study_col else ""
        if study and study in seen_studies:
            continue
        maybe_append(row)
        if len(specs) >= max_studies_per_locus:
            break

    if len(specs) >= max_studies_per_locus:
        return specs

    for _, row in eligible.iterrows():
        maybe_append(row)
        if len(specs) >= max_studies_per_locus:
            break
    return specs


def try_api_request(
    session: requests.Session,
    *,
    url: str,
    params: dict[str, Any],
    rsid: str,
    timeout_seconds: float,
) -> tuple[str, pd.DataFrame, str]:
    response = session.get(url, params=params, timeout=timeout_seconds)
    detail = compact_detail(response.text)

    try:
        payload = response.json()
    except ValueError:
        payload = None

    if payload is not None:
        detail = compact_detail(extract_message(payload) or detail)

    if response.status_code >= 400:
        return "error", pd.DataFrame(), f"HTTP {response.status_code}: {detail or 'request failed'}"

    if payload is None:
        return "error", pd.DataFrame(), "Non-JSON response from eQTL API."

    hits = parse_eqtl_api_payload(payload, lead_rsid=rsid)
    if not hits.empty:
        return "hits", hits, detail or f"{len(hits):,} associations returned."
    if payload_has_association_shape(payload):
        return "no_hits", hits, detail or "No associations returned."
    return "inconclusive", pd.DataFrame(), detail or "Unexpected response shape from eQTL API."


def lookup_eqtl_by_rsid(
    session: requests.Session,
    *,
    rsid: str,
    api_base_url: str,
    timeout_seconds: float,
    study_specs: list[dict[str, str]],
    pvalue_threshold: float | None,
) -> tuple[pd.DataFrame, str, str, str]:
    failures: list[str] = []
    all_hits: list[pd.DataFrame] = []
    completed_datasets = 0
    no_hit_datasets = 0

    for spec in study_specs:
        dataset_id = str(spec.get("dataset_id", "") or "").strip()
        if not dataset_id:
            continue
        dataset_failed = True
        params = {"rsid": rsid, "size": 1000}
        if pvalue_threshold is not None:
            params["p_upper"] = pvalue_threshold
        for api_version in ("v3", "v2"):
            outcome, hits, detail = try_api_request(
                session,
                url=f"{api_base_url}/{api_version}/datasets/{dataset_id}/associations",
                params=params,
                rsid=rsid,
                timeout_seconds=timeout_seconds,
            )
            mode = f"rsid_api_dataset:{api_version}:{dataset_id}"
            if outcome == "hits":
                dataset_failed = False
                completed_datasets += 1
                tagged_hits = hits.drop_duplicates().copy()
                tagged_hits["eqtl_dataset_id"] = dataset_id
                tagged_hits["eqtl_api_version"] = api_version
                all_hits.append(tagged_hits)
                break
            if outcome == "no_hits":
                dataset_failed = False
                completed_datasets += 1
                no_hit_datasets += 1
                break
            failures.append(f"{mode}: {detail}")
            if detail.startswith("HTTP 404:"):
                continue
            break
        if dataset_failed:
            continue

    if all_hits:
        merged = pd.concat(all_hits, ignore_index=True).drop_duplicates()
        status = "complete" if not failures else "partial_complete"
        detail = (
            f"Aggregated {len(merged):,} associations across {int(merged['eqtl_dataset_id'].nunique())} dataset(s)"
        )
        if no_hit_datasets:
            detail += f"; {no_hit_datasets} dataset(s) returned no hits"
        if failures:
            detail += f"; {len(failures)} dataset/version attempt(s) failed"
        return merged, status, "rsid_api_aggregate", compact_detail(detail)

    if completed_datasets:
        status = "no_hits" if not failures else "partial_complete"
        detail = f"No associations returned across {completed_datasets} completed dataset lookup(s)."
        if failures:
            detail += f" {len(failures)} dataset/version attempt(s) failed."
        return pd.DataFrame(), status, "rsid_api_aggregate", compact_detail(detail)

    raise RuntimeError(compact_detail("; ".join(failures) or "No successful rsID lookup attempt."))


def lookup_eqtl_by_variant(
    session: requests.Session,
    *,
    variant: str,
    api_base_url: str,
    timeout_seconds: float,
    study_specs: list[dict[str, str]],
    pvalue_threshold: float | None,
) -> tuple[pd.DataFrame, str, str, str]:
    failures: list[str] = []
    all_hits: list[pd.DataFrame] = []
    completed_datasets = 0
    no_hit_datasets = 0

    for spec in study_specs:
        dataset_id = str(spec.get("dataset_id", "") or "").strip()
        if not dataset_id:
            continue
        dataset_failed = True
        params = {"variant": variant, "size": 1000}
        if pvalue_threshold is not None:
            params["p_upper"] = pvalue_threshold
        for api_version in ("v3", "v2"):
            outcome, hits, detail = try_api_request(
                session,
                url=f"{api_base_url}/{api_version}/datasets/{dataset_id}/associations",
                params=params,
                rsid="",
                timeout_seconds=timeout_seconds,
            )
            mode = f"variant_api_dataset:{api_version}:{dataset_id}"
            if outcome == "hits":
                dataset_failed = False
                completed_datasets += 1
                tagged_hits = hits.drop_duplicates().copy()
                tagged_hits["eqtl_dataset_id"] = dataset_id
                tagged_hits["eqtl_api_version"] = api_version
                all_hits.append(tagged_hits)
                break
            if outcome == "no_hits":
                dataset_failed = False
                completed_datasets += 1
                no_hit_datasets += 1
                break
            failures.append(f"{mode}: {detail}")
            if detail.startswith("HTTP 404:"):
                continue
            break
        if dataset_failed:
            continue

    if all_hits:
        merged = pd.concat(all_hits, ignore_index=True).drop_duplicates()
        status = "complete" if not failures else "partial_complete"
        detail = (
            f"Aggregated {len(merged):,} associations across {int(merged['eqtl_dataset_id'].nunique())} dataset(s)"
        )
        if no_hit_datasets:
            detail += f"; {no_hit_datasets} dataset(s) returned no hits"
        if failures:
            detail += f"; {len(failures)} dataset/version attempt(s) failed"
        return merged, status, "variant_api_aggregate", compact_detail(detail)

    if completed_datasets:
        status = "no_hits" if not failures else "partial_complete"
        detail = f"No associations returned across {completed_datasets} completed dataset lookup(s)."
        if failures:
            detail += f" {len(failures)} dataset/version attempt(s) failed."
        return pd.DataFrame(), status, "variant_api_aggregate", compact_detail(detail)

    raise RuntimeError(compact_detail("; ".join(failures) or "No successful variant lookup attempt."))


def resolve_missing_rsids_via_variant_recoder(
    session: requests.Session,
    loci: pd.DataFrame,
    *,
    source_build: str | None,
    variant_recoder_url: str,
    timeout_seconds: float,
    batch_size: int = 50,
) -> tuple[dict[str, str], str | None]:
    if loci.empty:
        return {}, None
    if normalize_build_name(source_build) != "GRCh37":
        return {}, None

    pending_by_input: dict[str, list[str]] = {}
    for _, locus in loci.iterrows():
        if normalize_rsid_from_locus(locus):
            continue
        locus_id = str(locus.get("locus_id", "") or "").strip()
        if not locus_id:
            continue
        recoder_input = variant_recoder_input_from_locus(locus)
        if not recoder_input:
            continue
        pending_by_input.setdefault(recoder_input, []).append(locus_id)

    if not pending_by_input:
        return {}, None

    resolved_by_locus: dict[str, str] = {}
    pending_inputs = list(pending_by_input)
    try:
        for start in range(0, len(pending_inputs), batch_size):
            batch = pending_inputs[start : start + batch_size]
            response = session.post(
                variant_recoder_url,
                json={"ids": batch, "fields": "id"},
                timeout=timeout_seconds,
            )
            detail = compact_detail(response.text)
            if response.status_code >= 400:
                raise RuntimeError(f"HTTP {response.status_code}: {detail or 'variant_recoder request failed'}")
            try:
                payload = response.json()
            except ValueError as exc:
                raise RuntimeError("Non-JSON response from Ensembl variant_recoder.") from exc
            for recoder_input, rsid in extract_variant_recoder_rsids(payload).items():
                for locus_id in pending_by_input.get(recoder_input, []):
                    resolved_by_locus[locus_id] = rsid
    except Exception as exc:
        return {}, compact_detail(str(exc))

    return resolved_by_locus, None


def bridge_loci_to_eqtl_variants(
    session: requests.Session,
    loci: pd.DataFrame,
    *,
    source_build: str | None,
    target_build: str | None,
    grch37_rest_base_url: str,
    grch38_sequence_base_url: str,
    timeout_seconds: float,
) -> tuple[dict[str, str], dict[str, str]]:
    bridged: dict[str, str] = {}
    errors: dict[str, str] = {}

    source = normalize_build_name(source_build)
    target = normalize_build_name(target_build)
    if source != "GRCh37" or target != "GRCh38":
        return bridged, errors

    for _, locus in loci.iterrows():
        locus_id = str(locus.get("locus_id", "") or "").strip()
        if not locus_id:
            continue
        ref = str(locus.get("lead_ref", "") or "").strip()
        alt = str(locus.get("lead_alt", "") or "").strip()
        pos = pd.to_numeric(locus.get("lead_pos"), errors="coerce")
        if pd.isna(pos) or not ref or not alt:
            continue

        start = int(pos)
        end = start + len(ref) - 1
        region = f"X:{start}..{end}:1"
        try:
            map_response = session.get(
                f"{grch37_rest_base_url.rstrip('/')}/map/human/GRCh37/{region}/GRCh38",
                timeout=timeout_seconds,
            )
            if map_response.status_code >= 400:
                errors[locus_id] = f"Map request failed with HTTP {map_response.status_code}."
                continue
            mappings = map_response.json().get("mappings") or []
            if len(mappings) != 1:
                errors[locus_id] = "Map request did not return a single unique mapping."
                continue
            mapped = mappings[0].get("mapped") or {}
            mapped_start = pd.to_numeric(mapped.get("start"), errors="coerce")
            mapped_end = pd.to_numeric(mapped.get("end"), errors="coerce")
            mapped_chrom = str(mapped.get("seq_region_name", "") or "").strip()
            mapped_strand = int(mapped.get("strand") or 0)
            if pd.isna(mapped_start) or pd.isna(mapped_end) or not mapped_chrom:
                errors[locus_id] = "Mapped region was incomplete."
                continue
            if mapped_strand != 1:
                errors[locus_id] = "Mapped region was not on the forward strand."
                continue
            if int(mapped_end) - int(mapped_start) + 1 != len(ref):
                errors[locus_id] = "Mapped span length did not match the source reference length."
                continue

            sequence_response = session.get(
                f"{grch38_sequence_base_url.rstrip('/')}/{mapped_chrom}:{int(mapped_start)}..{int(mapped_end)}:1",
                timeout=timeout_seconds,
            )
            if sequence_response.status_code >= 400:
                errors[locus_id] = f"GRCh38 sequence fetch failed with HTTP {sequence_response.status_code}."
                continue
            mapped_ref = sequence_response.text.strip()
            if mapped_ref != ref:
                errors[locus_id] = "Mapped GRCh38 reference sequence did not match the source reference allele."
                continue

            variant = eqtl_variant_query_from_components(mapped_chrom, int(mapped_start), ref, alt)
            if variant:
                bridged[locus_id] = variant
        except Exception as exc:
            errors[locus_id] = compact_detail(str(exc))

    return bridged, errors


def main() -> None:
    parser = argparse.ArgumentParser(description="Query eQTL Catalogue for lead chrX loci with build-safe defaults.")
    parser.add_argument("--loci", default=str(PROCESSED_DIR / "x_loci.tsv.gz"))
    parser.add_argument("--eqtl-index", default=str(INTERIM_DIR / "eqtl_chrX_index.tsv"))
    parser.add_argument("--query-mode", choices=["rsid", "region_provisional"], default="rsid")
    parser.add_argument("--lead-build", default="GRCh37")
    parser.add_argument("--eqtl-build", default="GRCh38")
    parser.add_argument("--allow-unsafe-region-build", action="store_true")
    parser.add_argument("--api-base-url", default=EQTL_API_BASE_URL)
    parser.add_argument("--variant-recoder-url", default=ENSEMBL_GRCH37_VARIANT_RECODER_URL)
    parser.add_argument("--grch37-rest-base-url", default=ENSEMBL_GRCH37_REST_BASE_URL)
    parser.add_argument("--grch38-sequence-base-url", default=ENSEMBL_GRCH38_SEQUENCE_BASE_URL)
    parser.add_argument("--request-timeout-seconds", type=float, default=30.0)
    parser.add_argument("--variant-recoder-timeout-seconds", type=float, default=30.0)
    parser.add_argument("--variant-bridge-timeout-seconds", type=float, default=30.0)
    parser.add_argument(
        "--attempt-exact-variant-bridge",
        action="store_true",
        help="Experimental: attempt a GRCh37->GRCh38 exact-variant bridge for loci that still lack an rsID after variant_recoder recovery.",
    )
    parser.add_argument("--window-bp", type=int, default=1_000_000)
    parser.add_argument("--max-studies-per-locus", type=int, default=5)
    parser.add_argument(
        "--eqtl-pvalue-threshold",
        type=float,
        default=1e-5,
        help="Conservative raw p-value cutoff applied to API queries when per-dataset FDR calls are unavailable.",
    )
    parser.add_argument("--sleep-seconds", type=float, default=2.0)
    parser.add_argument("--out", default=str(INTERIM_DIR / "eqtl_hits.tsv.gz"))
    parser.add_argument("--summary-out", default=str(INTERIM_DIR / "eqtl_lookup_summary.tsv"))
    args = parser.parse_args()

    loci = pd.read_csv(args.loci, sep="\t", compression="gzip")
    eqtl_index = pd.read_csv(args.eqtl_index, sep="\t")

    if "chrX_supported" in eqtl_index.columns:
        eqtl_index = eqtl_index.loc[eqtl_index["chrX_supported"] == True].copy()
    if "priority" in eqtl_index.columns:
        eqtl_index = eqtl_index.sort_values(["priority", "study_id"], na_position="last")

    all_hits = []
    status_rows = []
    study_specs = choose_eqtl_study_specs(eqtl_index, args.max_studies_per_locus)

    path_col = None
    if args.query_mode == "region_provisional":
        if not shutil_which("tabix"):
            raise RuntimeError("Region-provisional mode requires 'tabix' on PATH.")
        path_col = choose_column(list(eqtl_index.columns), exact=["ftp_path_inferred"], contains_any=["ftp", "path", "url"])
        if path_col is None:
            raise RuntimeError("Could not find a path column in the eQTL index.")

    session = requests.Session()
    session.headers.update({"accept": "application/json", "content-type": "application/json"})

    recovered_rsids: dict[str, str] = {}
    variant_recoder_error = None
    bridged_variants: dict[str, str] = {}
    bridged_variant_errors: dict[str, str] = {}
    if args.query_mode == "rsid":
        recovered_rsids, variant_recoder_error = resolve_missing_rsids_via_variant_recoder(
            session,
            loci,
            source_build=args.lead_build,
            variant_recoder_url=args.variant_recoder_url,
            timeout_seconds=args.variant_recoder_timeout_seconds,
        )
        if recovered_rsids:
            print(f"[eqtl] recovered {len(recovered_rsids):,} missing rsIDs via Ensembl variant_recoder")
        elif variant_recoder_error:
            print(f"[eqtl] Ensembl variant_recoder unavailable: {variant_recoder_error}")

        if args.attempt_exact_variant_bridge:
            remaining_without_rsid = loci.loc[
                [
                    not normalize_rsid_from_locus(row) and str(row.get("locus_id", "") or "").strip() not in recovered_rsids
                    for _, row in loci.iterrows()
                ]
            ].copy()
            bridged_variants, bridged_variant_errors = bridge_loci_to_eqtl_variants(
                session,
                remaining_without_rsid,
                source_build=args.lead_build,
                target_build=args.eqtl_build,
                grch37_rest_base_url=args.grch37_rest_base_url,
                grch38_sequence_base_url=args.grch38_sequence_base_url,
                timeout_seconds=args.variant_bridge_timeout_seconds,
            )
            if bridged_variants:
                print(f"[eqtl] bridged {len(bridged_variants):,} remaining loci to exact GRCh38 variant queries")

    for _, locus in loci.iterrows():
        locus_id = str(locus["locus_id"])
        trait_id = str(locus["trait_id"])
        rsid = normalize_rsid_from_locus(locus)
        recovered_via_variant_recoder = False
        bridged_variant = None
        if not rsid:
            rsid = recovered_rsids.get(locus_id)
            recovered_via_variant_recoder = bool(rsid)
        if not rsid:
            bridged_variant = bridged_variants.get(locus_id)
        status = None
        hits = pd.DataFrame()

        if args.query_mode == "rsid":
            if rsid:
                print(f"[eqtl] {locus_id} <- {rsid}")
                try:
                    hits, query_status, query_mode, query_detail = lookup_eqtl_by_rsid(
                        session,
                        rsid=rsid,
                        api_base_url=args.api_base_url.rstrip("/"),
                        timeout_seconds=args.request_timeout_seconds,
                        study_specs=study_specs,
                        pvalue_threshold=args.eqtl_pvalue_threshold,
                    )
                    if recovered_via_variant_recoder:
                        query_mode = f"variant_recoder+{query_mode}"
                        query_detail = f"Resolved missing rsID via Ensembl GRCh37 variant_recoder. {query_detail}"
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode=query_mode,
                        query_status=query_status,
                        query_detail=query_detail,
                        lead_rsid=rsid,
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )
                except Exception as exc:
                    failure_mode = "variant_recoder_rsid_api" if recovered_via_variant_recoder else "rsid_api"
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode=failure_mode,
                        query_status="lookup_failed",
                        query_detail=compact_detail(str(exc)),
                        lead_rsid=rsid,
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )
            elif bridged_variant:
                print(f"[eqtl] {locus_id} <- {bridged_variant}")
                try:
                    hits, query_status, query_mode, query_detail = lookup_eqtl_by_variant(
                        session,
                        variant=bridged_variant,
                        api_base_url=args.api_base_url.rstrip("/"),
                        timeout_seconds=args.request_timeout_seconds,
                        study_specs=study_specs,
                        pvalue_threshold=args.eqtl_pvalue_threshold,
                    )
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode=f"variant_bridge+{query_mode}",
                        query_status=query_status,
                        query_detail=f"Bridged GRCh37 lead variant to a validated GRCh38 exact variant query. {query_detail}",
                        lead_rsid=locus.get("lead_rsid"),
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )
                except Exception as exc:
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode="variant_bridge+variant_api",
                        query_status="lookup_failed",
                        query_detail=compact_detail(str(exc)),
                        lead_rsid=locus.get("lead_rsid"),
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )
            else:
                bridge_detail = bridged_variant_errors.get(locus_id, "")
                if variant_recoder_error:
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode="variant_recoder_rsid_api",
                        query_status="lookup_failed",
                        query_detail=f"Lead variant has no usable rsID and Ensembl variant_recoder failed: {variant_recoder_error}",
                        lead_rsid=locus.get("lead_rsid"),
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )
                else:
                    query_mode = "variant_recoder_rsid_api"
                    detail = "Lead variant has no usable rsID after direct and Ensembl GRCh37 variant_recoder recovery."
                    if args.attempt_exact_variant_bridge:
                        query_mode = "variant_bridge+variant_recoder_rsid_api"
                        detail = f"{detail} No validated GRCh38 exact variant bridge was available."
                        if bridge_detail:
                            detail = f"{detail} {bridge_detail}"
                    else:
                        detail = f"{detail} Exact-variant build bridging is experimental and disabled by default."
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode=query_mode,
                        query_status="skipped_no_rsid",
                        query_detail=detail,
                        lead_rsid=locus.get("lead_rsid"),
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )
        else:
            build_safe = is_region_query_build_safe(args.lead_build, args.eqtl_build)
            if not build_safe and not args.allow_unsafe_region_build:
                status = build_locus_query_status(
                    locus_id=locus_id,
                    query_mode="region_provisional",
                    query_status="blocked_build_mismatch",
                    query_detail="Region query blocked because lead and eQTL builds do not match.",
                    lead_rsid=locus.get("lead_rsid"),
                    source_build=args.lead_build,
                    target_build=args.eqtl_build,
                )
            else:
                lead_pos = int(locus["lead_pos"])
                region = f"X:{max(1, lead_pos - args.window_bp)}-{lead_pos + args.window_bp}"
                chosen_studies = eqtl_index.head(args.max_studies_per_locus)
                failures = []

                for _, study in chosen_studies.iterrows():
                    source = study[path_col]
                    study_id = study.get("study_id") or study.get("study")
                    print(f"[eqtl] {locus_id} <- {study_id} {region}")
                    try:
                        region_hits = run_tabix_region(source, region)
                    except Exception as exc:
                        failures.append(f"{study_id}: {compact_detail(str(exc))}")
                        time.sleep(args.sleep_seconds)
                        continue

                    if region_hits.empty:
                        time.sleep(args.sleep_seconds)
                        continue

                    region_hits["study_id"] = study_id
                    all_hits.append(
                        region_hits.assign(
                            locus_id=locus_id,
                            trait_id=trait_id,
                            lead_rsid=rsid,
                            eqtl_query_mode="region_provisional",
                            eqtl_query_status="complete",
                        )
                    )
                    hits = pd.concat([hits, region_hits], ignore_index=True)
                    time.sleep(args.sleep_seconds)

                if not hits.empty:
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode="region_provisional",
                        query_status="complete",
                        query_detail="Provisional region query returned one or more hits.",
                        lead_rsid=rsid,
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )
                elif failures:
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode="region_provisional",
                        query_status="lookup_failed",
                        query_detail=compact_detail("; ".join(failures)),
                        lead_rsid=rsid,
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )
                else:
                    status = build_locus_query_status(
                        locus_id=locus_id,
                        query_mode="region_provisional",
                        query_status="no_hits",
                        query_detail="Provisional region query completed with no returned hits.",
                        lead_rsid=rsid,
                        source_build=args.lead_build,
                        target_build=args.eqtl_build,
                    )

        if args.query_mode == "rsid" and not hits.empty:
            hits = hits.assign(
                locus_id=locus_id,
                trait_id=trait_id,
                lead_rsid=rsid,
                eqtl_query_mode=status["query_mode"],
                eqtl_query_status=status["query_status"],
            )
            all_hits.append(hits)

        if status is None:
            status = build_locus_query_status(
                locus_id=locus_id,
                query_mode=args.query_mode,
                query_status="lookup_failed",
                query_detail="Lookup finished without a status record.",
                lead_rsid=rsid,
                source_build=args.lead_build,
                target_build=args.eqtl_build,
            )

        status["trait_id"] = trait_id
        status["n_hits"] = int(len(hits))
        status_rows.append(status)
        time.sleep(args.sleep_seconds)

    merged = pd.concat(all_hits, ignore_index=True) if all_hits else pd.DataFrame()
    summary = pd.DataFrame(status_rows)
    write_tsv(merged, Path(args.out))
    write_tsv(summary, Path(args.summary_out))
    print(f"[done] wrote {len(merged):,} rows to {args.out}")
    print(f"[done] wrote {len(summary):,} lookup status rows to {args.summary_out}")


def shutil_which(name: str) -> str | None:
    import shutil

    return shutil.which(name)


if __name__ == "__main__":
    main()
