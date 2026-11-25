from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from .io import coerce_bool, ensure_dir
from .scoring import effective_hq_case_count, stratum_label

SCHEMA_VERSION = 3
KNOWN_STRATA = ("both_sexes", "females", "males")


@dataclass(frozen=True)
class PanelStratumSpec:
    stratum: str
    label: str
    release_dir: Path


@dataclass(frozen=True)
class PanelBundleSpec:
    panel_id: str
    slug: str
    title: str
    description: str
    release_dir: Path | None = None
    featured: bool = True
    default_stratum: str = "both_sexes"
    strata: tuple[PanelStratumSpec, ...] = ()


@dataclass(frozen=True)
class PanelRelease:
    spec: PanelBundleSpec
    stratum_spec: PanelStratumSpec
    trait_scores: pd.DataFrame
    loci: pd.DataFrame
    gene_candidates: pd.DataFrame


def _detected_panel_strata(release_root: Path, panel_id: str) -> tuple[PanelStratumSpec, ...]:
    stratified_root = release_root / "stratified" / panel_id
    if not stratified_root.exists():
        return ()

    strata = []
    for stratum in KNOWN_STRATA:
        release_dir = stratified_root / stratum
        if (release_dir / "trait_scores.tsv").exists():
            strata.append(PanelStratumSpec(stratum=stratum, label=stratum_label(stratum), release_dir=release_dir))
    return tuple(strata)


def _panel_bundle_spec(
    release_root: Path,
    *,
    panel_id: str,
    slug: str,
    title: str,
    description: str,
    featured: bool,
) -> PanelBundleSpec:
    strata = _detected_panel_strata(release_root, panel_id)
    default_stratum = "both_sexes" if any(spec.stratum == "both_sexes" for spec in strata) else (strata[0].stratum if strata else "both_sexes")
    release_dir = None if strata else release_root / panel_id
    return PanelBundleSpec(
        panel_id=panel_id,
        slug=slug,
        title=title,
        description=description,
        release_dir=release_dir,
        featured=featured,
        default_stratum=default_stratum,
        strata=strata,
    )


def default_panel_specs(release_root: str | Path) -> list[PanelBundleSpec]:
    release_root = Path(release_root)
    return [
        _panel_bundle_spec(
            release_root,
            panel_id="panel_expanded",
            slug="broad-atlas",
            title="Broad Atlas",
            description="Curated all-purpose panel with broad trait coverage and low zero-locus drift.",
            featured=True,
        ),
        _panel_bundle_spec(
            release_root,
            panel_id="panel_mind_risk_core",
            slug="mind-and-behavior",
            title="Mind and Behavior",
            description="Focused behavior and cognition companion panel built around the strongest mind-adjacent chrX traits.",
            featured=True,
        ),
        _panel_bundle_spec(
            release_root,
            panel_id="panel_blood_biochemistry_core",
            slug="biochemistry-deep-dive",
            title="Biochemistry Deep Dive",
            description="High-yield blood biochemistry panel carved out of the independent-set discovery pool.",
            featured=True,
        ),
        _panel_bundle_spec(
            release_root,
            panel_id="panel_independent_set",
            slug="discovery-pool",
            title="Discovery Pool",
            description="Full Pan-UKB max independent-set build kept as a broad curation source rather than a featured panel.",
            featured=False,
        ),
    ]


def slugify(value: Any) -> str:
    text = "" if value is None else str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", "-", text)
    return text.strip("-") or "item"


def split_top_candidate_genes(value: Any) -> list[str]:
    if value is None or pd.isna(value):
        return []
    parts = [part.strip() for part in str(value).split(",")]
    return [part for part in parts if part]


def now_iso8601() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def normalize_json_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, dict):
        return {str(key): normalize_json_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize_json_value(item) for item in value]
    if isinstance(value, tuple):
        return [normalize_json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float, str)):
        if isinstance(value, float):
            if pd.isna(value):
                return None
            if value.is_integer():
                return int(value)
        return value
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        return normalize_json_value(value.item())
    return str(value)


def write_json(payload: dict[str, Any] | list[Any], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    normalized = normalize_json_value(payload)
    path.write_text(json.dumps(normalized, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def _bool_column(series: pd.Series) -> pd.Series:
    return series.map(lambda value: coerce_bool(value) is True)


def _int_column(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0).astype(int)


def _float_column(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _series_for_frame(df: pd.DataFrame, column: str, default: Any) -> pd.Series:
    if column in df.columns:
        return df[column]
    return pd.Series([default] * len(df), index=df.index)


def _sort_frame(df: pd.DataFrame, columns: list[str], ascending: list[bool]) -> pd.DataFrame:
    usable = [column for column in columns if column in df.columns]
    if not usable:
        return df
    usable_ascending = [ascending[columns.index(column)] for column in usable]
    return df.sort_values(usable, ascending=usable_ascending, na_position="last")


def panel_strata(spec: PanelBundleSpec) -> tuple[PanelStratumSpec, ...]:
    if spec.strata:
        return spec.strata
    if spec.release_dir is None:
        return ()
    return (
        PanelStratumSpec(
            stratum=spec.default_stratum,
            label=stratum_label(spec.default_stratum),
            release_dir=Path(spec.release_dir),
        ),
    )


def _default_panel_stratum(spec: PanelBundleSpec) -> PanelStratumSpec:
    strata = panel_strata(spec)
    if not strata:
        raise ValueError(f"Panel {spec.panel_id} has no release directories configured.")
    for stratum_spec in strata:
        if stratum_spec.stratum == spec.default_stratum:
            return stratum_spec
    return strata[0]


def _resolve_panel_stratum(spec: PanelBundleSpec, stratum: str | PanelStratumSpec | None) -> PanelStratumSpec:
    if isinstance(stratum, PanelStratumSpec):
        return stratum
    if stratum is None:
        return _default_panel_stratum(spec)
    for stratum_spec in panel_strata(spec):
        if stratum_spec.stratum == stratum:
            return stratum_spec
    raise KeyError(f"Unknown stratum for {spec.panel_id}: {stratum}")


def load_panel_release(spec: PanelBundleSpec, stratum: str | PanelStratumSpec | None = None) -> PanelRelease:
    stratum_spec = _resolve_panel_stratum(spec, stratum)

    trait_scores = pd.read_csv(stratum_spec.release_dir / "trait_scores.tsv", sep="\t").copy()
    loci = pd.read_csv(stratum_spec.release_dir / "x_loci.tsv", sep="\t").copy()
    gene_candidates = pd.read_csv(stratum_spec.release_dir / "x_gene_candidates.tsv", sep="\t").copy()

    trait_lookup_source = _series_for_frame(trait_scores, "eqtl_lookup_hit", False)
    trait_scores["eqtl_lookup_hit_bool"] = _bool_column(trait_lookup_source)
    trait_scores["eqtl_supported_bool"] = _bool_column(_series_for_frame(trait_scores, "eqtl_supported", False))
    trait_scores["n_loci"] = _int_column(_series_for_frame(trait_scores, "n_loci", 0))
    trait_scores["eqtl_lookup_hit_locus_count"] = _int_column(_series_for_frame(trait_scores, "eqtl_lookup_hit_locus_count", 0))
    trait_scores["eqtl_lookup_hit_count"] = _int_column(
        _series_for_frame(
            trait_scores,
            "eqtl_lookup_hit_count",
            _series_for_frame(trait_scores, "eqtl_total_hit_count", 0),
        )
    )
    trait_scores["eqtl_supported_locus_count"] = _int_column(_series_for_frame(trait_scores, "eqtl_supported_locus_count", 0))
    trait_scores["eqtl_total_hit_count"] = _int_column(_series_for_frame(trait_scores, "eqtl_total_hit_count", 0))
    trait_scores["x_evidence_score"] = _float_column(_series_for_frame(trait_scores, "x_evidence_score", 0)).fillna(0.0)
    view_stratum_source = trait_scores["view_stratum"] if "view_stratum" in trait_scores.columns else _series_for_frame(
        trait_scores,
        "pheno_sex",
        stratum_spec.stratum,
    )
    trait_scores["view_stratum"] = view_stratum_source.fillna(stratum_spec.stratum)
    stratum_label_source = (
        trait_scores["stratum_label"] if "stratum_label" in trait_scores.columns else trait_scores["view_stratum"].map(stratum_label)
    )
    trait_scores["stratum_label"] = stratum_label_source.fillna(trait_scores["view_stratum"].map(stratum_label))
    effective_hq_source = (
        trait_scores["effective_hq_case_count"]
        if "effective_hq_case_count" in trait_scores.columns
        else trait_scores.apply(effective_hq_case_count, axis=1)
    )
    trait_scores["effective_hq_case_count"] = _float_column(effective_hq_source)
    trait_scores["selection_used_fallback"] = _bool_column(_series_for_frame(trait_scores, "selection_used_fallback", False))
    trait_scores["slug"] = _series_for_frame(
        trait_scores,
        "query_id",
        _series_for_frame(trait_scores, "description", ""),
    ).map(slugify)

    loci_lookup_source = _series_for_frame(loci, "eqtl_lookup_hit", False)
    loci["eqtl_lookup_hit_bool"] = _bool_column(loci_lookup_source)
    loci["eqtl_supported_bool"] = _bool_column(_series_for_frame(loci, "eqtl_supported", False))
    loci["lead_neglog10_pvalue"] = _float_column(_series_for_frame(loci, "lead_neglog10_pvalue", 0))
    loci["eqtl_lookup_hit_gene_count"] = _int_column(_series_for_frame(loci, "eqtl_lookup_hit_gene_count", 0))
    loci["eqtl_supported_gene_count"] = _int_column(_series_for_frame(loci, "eqtl_supported_gene_count", 0))
    loci["eqtl_lookup_n_hits"] = _int_column(_series_for_frame(loci, "eqtl_lookup_n_hits", 0))

    gene_lookup_source = _series_for_frame(gene_candidates, "eqtl_lookup_hit", False)
    gene_candidates["eqtl_lookup_hit_bool"] = _bool_column(gene_lookup_source)
    gene_candidates["eqtl_supported_bool"] = _bool_column(_series_for_frame(gene_candidates, "eqtl_supported", False))
    gene_candidates["candidate_rank"] = _int_column(_series_for_frame(gene_candidates, "candidate_rank", 0))
    gene_candidates["distance_to_lead_bp"] = _int_column(_series_for_frame(gene_candidates, "distance_to_lead_bp", 0))
    gene_candidates["eqtl_lookup_hit_count"] = _int_column(_series_for_frame(gene_candidates, "eqtl_lookup_hit_count", 0))
    gene_candidates["eqtl_dataset_count"] = _int_column(_series_for_frame(gene_candidates, "eqtl_dataset_count", 0))
    gene_candidates["eqtl_gene_role"] = _series_for_frame(gene_candidates, "eqtl_gene_role", "candidate").fillna("candidate")

    return PanelRelease(spec=spec, stratum_spec=stratum_spec, trait_scores=trait_scores, loci=loci, gene_candidates=gene_candidates)


def _sort_traits(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(
        ["x_evidence_score", "n_loci", "eqtl_total_hit_count", "query_id"],
        ascending=[False, False, False, True],
        na_position="last",
    )


def _compact_trait_payload(row: pd.Series) -> dict[str, Any]:
    return {
        "trait_id": row.get("trait_id"),
        "query_id": row.get("query_id"),
        "slug": row.get("slug") or slugify(row.get("query_id") or row.get("description")),
        "description": row.get("description"),
        "coding_description": row.get("coding_description"),
        "domain": row.get("domain"),
        "view_stratum": row.get("view_stratum"),
        "stratum_label": row.get("stratum_label"),
        "effective_hq_case_count": row.get("effective_hq_case_count"),
        "selection_used_fallback": bool(row.get("selection_used_fallback", False)),
        "x_evidence_score": row.get("x_evidence_score"),
        "coverage_grade": row.get("coverage_grade"),
        "n_loci": row.get("n_loci"),
        "eqtl_lookup_hit": bool(row.get("eqtl_lookup_hit_bool", False)),
        "eqtl_lookup_hit_locus_count": row.get("eqtl_lookup_hit_locus_count", 0),
        "eqtl_lookup_hit_count": row.get("eqtl_lookup_hit_count", row.get("eqtl_total_hit_count", 0)),
        "eqtl_supported": bool(row.get("eqtl_supported_bool", False)),
        "eqtl_supported_locus_count": row.get("eqtl_supported_locus_count", 0),
        "eqtl_total_hit_count": row.get("eqtl_total_hit_count", 0),
        "top_candidate_genes": split_top_candidate_genes(row.get("top_candidate_genes")),
        "confidence_notes": row.get("confidence_notes"),
        "eqtl_lookup_note": row.get("eqtl_lookup_note"),
    }


def build_traits_payload(panel: PanelRelease) -> list[dict[str, Any]]:
    ranked = _sort_traits(panel.trait_scores)
    return normalize_json_value([_compact_trait_payload(row) for _, row in ranked.iterrows()])


def _domain_summary(panel: PanelRelease) -> list[dict[str, Any]]:
    rows = []
    grouped = panel.trait_scores.groupby("domain", dropna=False)
    for domain, group in grouped:
        rows.append(
            {
                "domain": domain,
                "trait_count": int(len(group)),
                "lookup_hit_trait_count": int(group["eqtl_lookup_hit_bool"].sum()),
                "supported_trait_count": int(group["eqtl_supported_bool"].sum()),
                "zero_locus_trait_count": int(group["n_loci"].eq(0).sum()),
                "locus_count": int(group["n_loci"].sum()),
            }
        )
    rows.sort(key=lambda row: (-row["supported_trait_count"], -row["locus_count"], str(row["domain"])))
    return rows


def panel_metrics(panel: PanelRelease) -> dict[str, Any]:
    return {
        "trait_count": int(len(panel.trait_scores)),
        "lookup_hit_trait_count": int(panel.trait_scores["eqtl_lookup_hit_bool"].sum()),
        "supported_trait_count": int(panel.trait_scores["eqtl_supported_bool"].sum()),
        "zero_locus_trait_count": int(panel.trait_scores["n_loci"].eq(0).sum()),
        "locus_count": int(len(panel.loci)),
    }


def curated_trait_count(panel: PanelRelease) -> int:
    trait_scores = panel.trait_scores
    if trait_scores.empty:
        return 0
    exact = trait_scores.loc[
        trait_scores["view_stratum"].fillna("").eq(panel.stratum_spec.stratum)
        & ~trait_scores["selection_used_fallback"].fillna(False)
    ].copy()
    if exact.empty:
        return 0
    key = "query_id" if "query_id" in exact.columns else "trait_id"
    return int(exact[key].dropna().astype(str).nunique())


def _stratum_manifest_entry(panel: PanelRelease) -> dict[str, Any]:
    metrics = panel_metrics(panel)
    return {
        "label": panel.stratum_spec.label,
        "curated_trait_count": curated_trait_count(panel),
        **metrics,
        "files": {
            "summary": f"panels/{panel.spec.panel_id}/{panel.stratum_spec.stratum}/summary.json",
            "traits": f"panels/{panel.spec.panel_id}/{panel.stratum_spec.stratum}/traits.json",
            "trait_details_dir": f"panels/{panel.spec.panel_id}/{panel.stratum_spec.stratum}/traits",
        },
    }


def manifest_entry(panel: PanelRelease) -> dict[str, Any]:
    stratum_entry = _stratum_manifest_entry(panel)
    return {
        "panel_id": panel.spec.panel_id,
        "slug": panel.spec.slug,
        "title": panel.spec.title,
        "featured": panel.spec.featured,
        "description": panel.spec.description,
        "default_stratum": panel.stratum_spec.stratum,
        "available_strata": [panel.stratum_spec.stratum],
        **panel_metrics(panel),
        "files": stratum_entry["files"],
        "strata": {panel.stratum_spec.stratum: stratum_entry},
    }


def _manifest_panel_entry(
    spec: PanelBundleSpec,
    releases: list[PanelRelease],
    *,
    min_curated_traits_per_stratum: int = 1,
) -> dict[str, Any]:
    if not releases:
        raise ValueError(f"Panel {spec.panel_id} produced no releases.")

    by_stratum = {panel.stratum_spec.stratum: panel for panel in releases}
    exposed = {
        stratum: panel
        for stratum, panel in by_stratum.items()
        if curated_trait_count(panel) >= max(1, min_curated_traits_per_stratum)
    }
    if not exposed:
        exposed = {releases[0].stratum_spec.stratum: releases[0]}

    default_stratum = spec.default_stratum if spec.default_stratum in exposed else next(iter(exposed))
    if default_stratum not in exposed and "both_sexes" in exposed:
        default_stratum = "both_sexes"
    default_release = exposed[default_stratum]
    default_entry = _stratum_manifest_entry(default_release)

    return {
        "panel_id": spec.panel_id,
        "slug": spec.slug,
        "title": spec.title,
        "featured": spec.featured,
        "description": spec.description,
        "default_stratum": default_stratum,
        "available_strata": list(exposed.keys()),
        **panel_metrics(default_release),
        "files": default_entry["files"],
        "strata": {stratum: _stratum_manifest_entry(panel) for stratum, panel in exposed.items()},
    }


def build_panel_summary_payload(panel: PanelRelease, *, generated_at: str) -> dict[str, Any]:
    ranked = build_traits_payload(panel)
    unsupported = [row for row in ranked if not row["eqtl_supported"]]
    return normalize_json_value(
        {
            "schema_version": SCHEMA_VERSION,
            "generated_at": generated_at,
            "view": {
                "panel_id": panel.spec.panel_id,
                "stratum": panel.stratum_spec.stratum,
                "label": panel.stratum_spec.label,
            },
            "panel": {
                "panel_id": panel.spec.panel_id,
                "slug": panel.spec.slug,
                "title": panel.spec.title,
                "featured": panel.spec.featured,
                "description": panel.spec.description,
            },
            "metrics": panel_metrics(panel),
            "domain_summary": _domain_summary(panel),
            "top_traits": ranked[:20],
            "unsupported_traits": unsupported[:20],
        }
    )


def _trait_metadata(row: pd.Series) -> dict[str, Any]:
    extra_fields = [
        "trait_type",
        "phenocode",
        "pheno_sex",
        "view_stratum",
        "stratum_label",
        "coding",
        "modifier",
        "category",
        "description_more",
        "coding_description",
        "num_pops_pass_qc",
        "in_max_independent_set",
        "tier",
        "selection_rank",
        "selection_score",
        "selection_used_fallback",
        "requested_pheno_sex",
        "lead_neglog10_pvalue",
        "any_nonpar_locus",
        "par_only_signal",
        "effective_hq_case_count",
        "n_cases_hq_cohort_both_sexes",
        "n_cases_hq_cohort_females",
        "n_cases_hq_cohort_males",
    ]
    payload = _compact_trait_payload(row)
    for field in extra_fields:
        if field in row.index:
            payload[field] = row.get(field)
    return payload


def _candidate_gene_payload(row: pd.Series) -> dict[str, Any]:
    return {
        "candidate_rank": row.get("candidate_rank"),
        "gene_id": row.get("gene_id"),
        "gene_name": row.get("gene_name"),
        "gene_biotype": row.get("gene_biotype"),
        "mapping_relation": row.get("mapping_relation"),
        "distance_to_lead_bp": row.get("distance_to_lead_bp"),
        "eqtl_gene_role": row.get("eqtl_gene_role"),
        "eqtl_lookup_hit": bool(row.get("eqtl_lookup_hit_bool", False)),
        "eqtl_lookup_hit_count": row.get("eqtl_lookup_hit_count"),
        "eqtl_supported": bool(row.get("eqtl_supported_bool", False)),
        "eqtl_study_count": row.get("eqtl_study_count"),
        "eqtl_dataset_count": row.get("eqtl_dataset_count"),
        "best_eqtl_pvalue": row.get("best_eqtl_pvalue"),
    }


def build_trait_detail_payload(
    panel: PanelRelease,
    trait_id: str,
    *,
    generated_at: str,
) -> dict[str, Any]:
    trait_rows = panel.trait_scores.loc[panel.trait_scores["trait_id"] == trait_id]
    if trait_rows.empty:
        raise KeyError(f"Unknown trait_id: {trait_id}")
    trait_row = trait_rows.iloc[0]
    loci = panel.loci.loc[panel.loci["trait_id"] == trait_id].copy()
    loci = _sort_frame(loci, ["lead_neglog10_pvalue", "locus_start"], [False, True])

    loci_payload = []
    for _, locus_row in loci.iterrows():
        locus_id = locus_row["locus_id"]
        gene_rows = panel.gene_candidates.loc[panel.gene_candidates["locus_id"] == locus_id].copy()
        gene_rows = _sort_frame(gene_rows, ["candidate_rank", "distance_to_lead_bp", "gene_name"], [True, True, True])
        loci_payload.append(
            {
                "locus_id": locus_row.get("locus_id"),
                "x_region": locus_row.get("x_region"),
                "locus_start": locus_row.get("locus_start"),
                "locus_end": locus_row.get("locus_end"),
                "lead_pos": locus_row.get("lead_pos"),
                "lead_rsid": locus_row.get("lead_rsid"),
                "lead_varid": locus_row.get("lead_varid"),
                "lead_neglog10_pvalue": locus_row.get("lead_neglog10_pvalue"),
                "top_gene_id": locus_row.get("top_gene_id"),
                "top_gene_name": locus_row.get("top_gene_name"),
                "best_mapping_relation": locus_row.get("best_mapping_relation"),
                "eqtl_lookup_hit": bool(locus_row.get("eqtl_lookup_hit_bool", False)),
                "eqtl_lookup_hit_gene_count": locus_row.get("eqtl_lookup_hit_gene_count"),
                "eqtl_supported": bool(locus_row.get("eqtl_supported_bool", False)),
                "eqtl_supported_gene_count": locus_row.get("eqtl_supported_gene_count"),
                "eqtl_lookup_status": locus_row.get("eqtl_lookup_status"),
                "eqtl_lookup_mode": locus_row.get("eqtl_lookup_mode"),
                "eqtl_lookup_n_hits": locus_row.get("eqtl_lookup_n_hits"),
                "candidate_genes": [_candidate_gene_payload(gene_row) for _, gene_row in gene_rows.iterrows()],
            }
        )

    return normalize_json_value(
        {
            "schema_version": SCHEMA_VERSION,
            "generated_at": generated_at,
            "view": {
                "panel_id": panel.spec.panel_id,
                "stratum": panel.stratum_spec.stratum,
                "label": panel.stratum_spec.label,
            },
            "panel": {
                "panel_id": panel.spec.panel_id,
                "slug": panel.spec.slug,
                "title": panel.spec.title,
                "featured": panel.spec.featured,
            },
            "trait": _trait_metadata(trait_row),
            "loci": loci_payload,
        }
    )


def export_panel_bundle(panel: PanelRelease, outdir: str | Path, *, generated_at: str) -> dict[str, Any]:
    outdir = Path(outdir)
    panel_dir = ensure_dir(outdir / "panels" / panel.spec.panel_id / panel.stratum_spec.stratum)
    trait_dir = ensure_dir(panel_dir / "traits")

    summary_payload = build_panel_summary_payload(panel, generated_at=generated_at)
    traits_payload = build_traits_payload(panel)

    write_json(summary_payload, panel_dir / "summary.json")
    write_json(traits_payload, panel_dir / "traits.json")
    for trait_row in traits_payload:
        write_json(
            build_trait_detail_payload(panel, trait_row["trait_id"], generated_at=generated_at),
            trait_dir / f'{trait_row["trait_id"]}.json',
        )

    return _stratum_manifest_entry(panel)


def export_frontend_bundle(
    panel_specs: list[PanelBundleSpec],
    outdir: str | Path,
    *,
    generated_at: str | None = None,
    min_curated_traits_per_stratum: int = 1,
) -> dict[str, Any]:
    outdir = ensure_dir(outdir)
    generated_at = generated_at or now_iso8601()

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": generated_at,
        "panels": [],
    }
    for spec in panel_specs:
        releases = []
        for stratum_spec in panel_strata(spec):
            panel = load_panel_release(spec, stratum_spec)
            export_panel_bundle(panel, outdir, generated_at=generated_at)
            releases.append(panel)
        manifest["panels"].append(
            _manifest_panel_entry(
                spec,
                releases,
                min_curated_traits_per_stratum=min_curated_traits_per_stratum,
            )
        )

    write_json(manifest, Path(outdir) / "manifest.json")
    return manifest
