from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd

from .io import coerce_bool, ensure_dir

SCHEMA_VERSION = 1


@dataclass(frozen=True)
class PanelBundleSpec:
    panel_id: str
    slug: str
    title: str
    description: str
    release_dir: Path
    featured: bool = True


@dataclass(frozen=True)
class PanelRelease:
    spec: PanelBundleSpec
    trait_scores: pd.DataFrame
    loci: pd.DataFrame
    gene_candidates: pd.DataFrame


def default_panel_specs(release_root: str | Path) -> list[PanelBundleSpec]:
    release_root = Path(release_root)
    return [
        PanelBundleSpec(
            panel_id="panel_expanded",
            slug="broad-atlas",
            title="Broad Atlas",
            description="Curated all-purpose panel with broad trait coverage and low zero-locus drift.",
            release_dir=release_root / "panel_expanded",
            featured=True,
        ),
        PanelBundleSpec(
            panel_id="panel_mind_risk_core",
            slug="mind-and-risk",
            title="Mind And Risk",
            description="Focused behavior and cognition companion panel built around the strongest mind-adjacent chrX traits.",
            release_dir=release_root / "panel_mind_risk_core",
            featured=True,
        ),
        PanelBundleSpec(
            panel_id="panel_blood_biochemistry_core",
            slug="biochemistry-deep-dive",
            title="Biochemistry Deep Dive",
            description="High-yield blood biochemistry panel carved out of the independent-set discovery pool.",
            release_dir=release_root / "panel_blood_biochemistry_core",
            featured=True,
        ),
        PanelBundleSpec(
            panel_id="panel_independent_set",
            slug="discovery-pool",
            title="Discovery Pool",
            description="Full Pan-UKB max independent-set build kept as a broad curation source rather than a featured panel.",
            release_dir=release_root / "panel_independent_set",
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


def load_panel_release(spec: PanelBundleSpec) -> PanelRelease:
    trait_scores = pd.read_csv(spec.release_dir / "trait_scores.tsv", sep="\t").copy()
    loci = pd.read_csv(spec.release_dir / "x_loci.tsv", sep="\t").copy()
    gene_candidates = pd.read_csv(spec.release_dir / "x_gene_candidates.tsv", sep="\t").copy()

    trait_scores["eqtl_supported_bool"] = _bool_column(trait_scores["eqtl_supported"])
    trait_scores["n_loci"] = _int_column(trait_scores["n_loci"])
    trait_scores["eqtl_supported_locus_count"] = _int_column(trait_scores.get("eqtl_supported_locus_count", 0))
    trait_scores["eqtl_total_hit_count"] = _int_column(trait_scores.get("eqtl_total_hit_count", 0))
    trait_scores["x_evidence_score"] = _float_column(trait_scores.get("x_evidence_score", 0)).fillna(0.0)
    trait_scores["slug"] = trait_scores.get("query_id", trait_scores.get("description")).map(slugify)

    loci["eqtl_supported_bool"] = _bool_column(loci.get("eqtl_supported", False))
    loci["lead_neglog10_pvalue"] = _float_column(loci.get("lead_neglog10_pvalue"))
    loci["eqtl_lookup_n_hits"] = _int_column(loci.get("eqtl_lookup_n_hits", 0))

    gene_candidates["eqtl_supported_bool"] = _bool_column(gene_candidates.get("eqtl_supported", False))
    gene_candidates["candidate_rank"] = _int_column(gene_candidates.get("candidate_rank", 0))
    gene_candidates["distance_to_lead_bp"] = _int_column(gene_candidates.get("distance_to_lead_bp", 0))

    return PanelRelease(spec=spec, trait_scores=trait_scores, loci=loci, gene_candidates=gene_candidates)


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
        "domain": row.get("domain"),
        "x_evidence_score": row.get("x_evidence_score"),
        "coverage_grade": row.get("coverage_grade"),
        "n_loci": row.get("n_loci"),
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
        "supported_trait_count": int(panel.trait_scores["eqtl_supported_bool"].sum()),
        "zero_locus_trait_count": int(panel.trait_scores["n_loci"].eq(0).sum()),
        "locus_count": int(len(panel.loci)),
    }


def manifest_entry(panel: PanelRelease) -> dict[str, Any]:
    metrics = panel_metrics(panel)
    return {
        "panel_id": panel.spec.panel_id,
        "slug": panel.spec.slug,
        "title": panel.spec.title,
        "featured": panel.spec.featured,
        "description": panel.spec.description,
        **metrics,
        "files": {
            "summary": f"panels/{panel.spec.panel_id}/summary.json",
            "traits": f"panels/{panel.spec.panel_id}/traits.json",
            "trait_details_dir": f"panels/{panel.spec.panel_id}/traits",
        },
    }


def build_panel_summary_payload(panel: PanelRelease, *, generated_at: str) -> dict[str, Any]:
    ranked = build_traits_payload(panel)
    unsupported = [row for row in ranked if not row["eqtl_supported"]]
    return normalize_json_value({
        "schema_version": SCHEMA_VERSION,
        "generated_at": generated_at,
        "panel": manifest_entry(panel),
        "metrics": panel_metrics(panel),
        "domain_summary": _domain_summary(panel),
        "top_traits": ranked[:20],
        "unsupported_traits": unsupported[:20],
    })


def _trait_metadata(row: pd.Series) -> dict[str, Any]:
    extra_fields = [
        "trait_type",
        "phenocode",
        "pheno_sex",
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
        "lead_neglog10_pvalue",
        "any_nonpar_locus",
        "par_only_signal",
        "n_cases_hq_cohort_both_sexes",
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
        "eqtl_supported": bool(row.get("eqtl_supported_bool", False)),
        "eqtl_study_count": row.get("eqtl_study_count"),
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
    loci = loci.sort_values(["lead_neglog10_pvalue", "locus_start"], ascending=[False, True], na_position="last")

    loci_payload = []
    for _, locus_row in loci.iterrows():
        locus_id = locus_row["locus_id"]
        gene_rows = panel.gene_candidates.loc[panel.gene_candidates["locus_id"] == locus_id].copy()
        gene_rows = gene_rows.sort_values(["candidate_rank", "distance_to_lead_bp", "gene_name"], ascending=[True, True, True], na_position="last")
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
                "eqtl_supported": bool(locus_row.get("eqtl_supported_bool", False)),
                "eqtl_lookup_status": locus_row.get("eqtl_lookup_status"),
                "eqtl_lookup_mode": locus_row.get("eqtl_lookup_mode"),
                "eqtl_lookup_n_hits": locus_row.get("eqtl_lookup_n_hits"),
                "candidate_genes": [_candidate_gene_payload(gene_row) for _, gene_row in gene_rows.iterrows()],
            }
        )

    return normalize_json_value({
        "schema_version": SCHEMA_VERSION,
        "generated_at": generated_at,
        "panel": {
            "panel_id": panel.spec.panel_id,
            "slug": panel.spec.slug,
            "title": panel.spec.title,
            "featured": panel.spec.featured,
        },
        "trait": _trait_metadata(trait_row),
        "loci": loci_payload,
    })


def export_panel_bundle(panel: PanelRelease, outdir: str | Path, *, generated_at: str) -> dict[str, Any]:
    outdir = Path(outdir)
    panel_dir = ensure_dir(outdir / "panels" / panel.spec.panel_id)
    trait_dir = ensure_dir(panel_dir / "traits")

    summary_payload = build_panel_summary_payload(panel, generated_at=generated_at)
    traits_payload = build_traits_payload(panel)

    write_json(summary_payload, panel_dir / "summary.json")
    write_json(traits_payload, panel_dir / "traits.json")
    for trait_row in traits_payload:
        write_json(build_trait_detail_payload(panel, trait_row["trait_id"], generated_at=generated_at), trait_dir / f'{trait_row["trait_id"]}.json')

    return manifest_entry(panel)


def export_frontend_bundle(
    panel_specs: list[PanelBundleSpec],
    outdir: str | Path,
    *,
    generated_at: str | None = None,
) -> dict[str, Any]:
    outdir = ensure_dir(outdir)
    generated_at = generated_at or now_iso8601()

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": generated_at,
        "panels": [],
    }
    for spec in panel_specs:
        panel = load_panel_release(spec)
        manifest["panels"].append(export_panel_bundle(panel, outdir, generated_at=generated_at))

    write_json(manifest, Path(outdir) / "manifest.json")
    return manifest
