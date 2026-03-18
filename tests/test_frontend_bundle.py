import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

from xatlas.frontend_bundle import (
    PanelBundleSpec,
    build_panel_summary_payload,
    build_trait_detail_payload,
    build_traits_payload,
    export_frontend_bundle,
    load_panel_release,
    manifest_entry,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def _write_release_dir(base: Path, trait_rows: list[dict], loci_rows: list[dict], gene_rows: list[dict]) -> Path:
    base.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(trait_rows).to_csv(base / "trait_scores.tsv", sep="\t", index=False)
    pd.DataFrame(loci_rows).to_csv(base / "x_loci.tsv", sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(base / "x_gene_candidates.tsv", sep="\t", index=False)
    return base


def _fixture_release(tmp_path: Path, name: str, *, featured: bool = True) -> PanelBundleSpec:
    release_dir = _write_release_dir(
        tmp_path / name,
        [
            {
                "trait_id": "trait_a",
                "query_id": "trait_a",
                "description": "Trait A",
                "domain": "alpha",
                "x_evidence_score": 90.0,
                "coverage_grade": "A",
                "n_loci": 2,
                "eqtl_supported": "True",
                "eqtl_supported_locus_count": 1,
                "eqtl_total_hit_count": 7,
                "top_candidate_genes": "GENE1,GENE2",
                "confidence_notes": "Strong signal.",
                "eqtl_lookup_note": "Completed.",
                "trait_type": "continuous",
                "phenocode": "100",
                "pheno_sex": "both_sexes",
                "coding": "",
                "modifier": "",
                "category": "Alpha",
                "description_more": "More text",
                "coding_description": "",
                "num_pops_pass_qc": 3,
                "in_max_independent_set": True,
                "tier": "A",
                "selection_rank": 1,
                "selection_score": 50,
                "lead_neglog10_pvalue": 12.3,
                "any_nonpar_locus": True,
                "par_only_signal": False,
                "n_cases_hq_cohort_both_sexes": 100000,
            },
            {
                "trait_id": "trait_b",
                "query_id": "trait_b",
                "description": "Trait B",
                "domain": "beta",
                "x_evidence_score": 55.0,
                "coverage_grade": "B",
                "n_loci": 0,
                "eqtl_supported": "False",
                "eqtl_supported_locus_count": 0,
                "eqtl_total_hit_count": 0,
                "top_candidate_genes": "",
                "confidence_notes": "",
                "eqtl_lookup_note": "",
                "trait_type": "categorical",
                "phenocode": "200",
                "pheno_sex": "both_sexes",
                "coding": "",
                "modifier": "",
                "category": "Beta",
                "description_more": "",
                "coding_description": "",
                "num_pops_pass_qc": 2,
                "in_max_independent_set": False,
                "tier": "B",
                "selection_rank": 2,
                "selection_score": 20,
                "lead_neglog10_pvalue": 0,
                "any_nonpar_locus": False,
                "par_only_signal": False,
                "n_cases_hq_cohort_both_sexes": 5000,
            },
        ],
        [
            {
                "trait_id": "trait_a",
                "locus_id": "locus_2",
                "x_region": "nonPAR",
                "locus_start": 20,
                "locus_end": 30,
                "lead_pos": 25,
                "lead_rsid": "rs2",
                "lead_varid": "X_25_A_G",
                "lead_neglog10_pvalue": 11.0,
                "top_gene_id": "ENSG2",
                "top_gene_name": "GENE2",
                "best_mapping_relation": "nearby_window",
                "eqtl_supported": False,
                "eqtl_lookup_status": "no_hits",
                "eqtl_lookup_mode": "rsid_api",
                "eqtl_lookup_n_hits": 0,
            },
            {
                "trait_id": "trait_a",
                "locus_id": "locus_1",
                "x_region": "nonPAR",
                "locus_start": 10,
                "locus_end": 15,
                "lead_pos": 12,
                "lead_rsid": "rs1",
                "lead_varid": "X_12_A_G",
                "lead_neglog10_pvalue": 15.0,
                "top_gene_id": "ENSG1",
                "top_gene_name": "GENE1",
                "best_mapping_relation": "lead_overlap",
                "eqtl_supported": True,
                "eqtl_lookup_status": "complete",
                "eqtl_lookup_mode": "rsid_api_dataset:v3:QTD000116",
                "eqtl_lookup_n_hits": 3,
            },
        ],
        [
            {
                "trait_id": "trait_a",
                "locus_id": "locus_1",
                "candidate_rank": 1,
                "gene_id": "ENSG1",
                "gene_name": "GENE1",
                "gene_biotype": "protein_coding",
                "mapping_relation": "lead_overlap",
                "distance_to_lead_bp": 0,
                "eqtl_supported": True,
                "eqtl_study_count": 2,
                "best_eqtl_pvalue": 1e-8,
            },
            {
                "trait_id": "trait_a",
                "locus_id": "locus_1",
                "candidate_rank": 2,
                "gene_id": "ENSGX",
                "gene_name": "GENEX",
                "gene_biotype": "lncRNA",
                "mapping_relation": "nearby_window",
                "distance_to_lead_bp": 250,
                "eqtl_supported": False,
                "eqtl_study_count": 0,
                "best_eqtl_pvalue": None,
            },
            {
                "trait_id": "trait_a",
                "locus_id": "locus_2",
                "candidate_rank": 1,
                "gene_id": "ENSG2",
                "gene_name": "GENE2",
                "gene_biotype": "protein_coding",
                "mapping_relation": "nearby_window",
                "distance_to_lead_bp": 40,
                "eqtl_supported": False,
                "eqtl_study_count": 0,
                "best_eqtl_pvalue": None,
            },
        ],
    )
    return PanelBundleSpec(name, f"{name}-slug", name.title(), f"{name} description", release_dir, featured=featured)


def test_manifest_entry_respects_featured_flag(tmp_path: Path) -> None:
    featured_panel = load_panel_release(_fixture_release(tmp_path, "featured_panel", featured=True))
    discovery_panel = load_panel_release(_fixture_release(tmp_path, "discovery_panel", featured=False))

    featured_manifest = manifest_entry(featured_panel)
    discovery_manifest = manifest_entry(discovery_panel)

    assert featured_manifest["featured"] is True
    assert discovery_manifest["featured"] is False
    assert discovery_manifest["files"]["summary"] == "panels/discovery_panel/summary.json"


def test_traits_payload_normalizes_values_and_sorts_stably(tmp_path: Path) -> None:
    panel = load_panel_release(_fixture_release(tmp_path, "sort_panel"))
    traits_payload = build_traits_payload(panel)

    assert traits_payload[0]["trait_id"] == "trait_a"
    assert traits_payload[0]["eqtl_supported"] is True
    assert traits_payload[0]["top_candidate_genes"] == ["GENE1", "GENE2"]
    assert traits_payload[1]["eqtl_supported"] is False


def test_trait_detail_nests_loci_and_candidate_genes(tmp_path: Path) -> None:
    panel = load_panel_release(_fixture_release(tmp_path, "detail_panel"))
    payload = build_trait_detail_payload(panel, "trait_a", generated_at="2026-03-18T00:00:00Z")

    assert payload["trait"]["trait_id"] == "trait_a"
    assert [row["locus_id"] for row in payload["loci"]] == ["locus_1", "locus_2"]
    assert payload["loci"][0]["candidate_genes"][0]["gene_id"] == "ENSG1"
    assert payload["loci"][0]["candidate_genes"][1]["best_eqtl_pvalue"] is None


def test_panel_summary_reports_metrics_and_domain_summary(tmp_path: Path) -> None:
    panel = load_panel_release(_fixture_release(tmp_path, "summary_panel"))
    payload = build_panel_summary_payload(panel, generated_at="2026-03-18T00:00:00Z")

    assert payload["metrics"]["supported_trait_count"] == 1
    assert payload["metrics"]["zero_locus_trait_count"] == 1
    assert payload["domain_summary"][0]["trait_count"] == 1
    assert payload["top_traits"][0]["trait_id"] == "trait_a"
    assert payload["unsupported_traits"][0]["trait_id"] == "trait_b"


def test_export_frontend_bundle_and_script_smoke(tmp_path: Path) -> None:
    panel_specs = [
        _fixture_release(tmp_path, "panel_expanded", featured=True),
        _fixture_release(tmp_path, "panel_mind_risk_core", featured=True),
        _fixture_release(tmp_path, "panel_blood_biochemistry_core", featured=True),
        _fixture_release(tmp_path, "panel_independent_set", featured=False),
    ]
    outdir = tmp_path / "frontend"

    manifest = export_frontend_bundle(panel_specs, outdir, generated_at="2026-03-18T00:00:00Z")
    assert len(manifest["panels"]) == 4
    assert (outdir / "manifest.json").exists()
    assert (outdir / "panels" / "panel_expanded" / "summary.json").exists()
    assert (outdir / "panels" / "panel_expanded" / "traits.json").exists()
    assert (outdir / "panels" / "panel_expanded" / "traits" / "trait_a.json").exists()

    script = REPO_ROOT / "scripts" / "14_export_frontend_bundle.py"
    cmd = [
        sys.executable,
        str(script),
        "--outdir",
        str(tmp_path / "frontend_script"),
        "--expanded-release-dir",
        str(tmp_path / "panel_expanded"),
        "--mind-release-dir",
        str(tmp_path / "panel_mind_risk_core"),
        "--biochem-release-dir",
        str(tmp_path / "panel_blood_biochemistry_core"),
        "--independent-release-dir",
        str(tmp_path / "panel_independent_set"),
    ]
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)

    manifest_payload = json.loads((tmp_path / "frontend_script" / "manifest.json").read_text(encoding="utf-8"))
    assert len(manifest_payload["panels"]) == 4
    assert next(panel for panel in manifest_payload["panels"] if panel["panel_id"] == "panel_independent_set")["featured"] is False
