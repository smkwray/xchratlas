import runpy
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
EXPORT_MOD = runpy.run_path(str(REPO_ROOT / "scripts/08_export_release_tables.py"))


def test_build_trait_gene_rollup_adds_eqtl_summary_counts():
    build_trait_gene_rollup = EXPORT_MOD["build_trait_gene_rollup"]
    gene_candidates = pd.DataFrame(
        [
            {
                "trait_id": "trait1",
                "locus_id": "locus1",
                "candidate_rank": 1,
                "gene_id": "ENSG1",
                "gene_name": "GENE1",
                "eqtl_gene_role": "candidate",
                "eqtl_lookup_hit": True,
                "mapping_relation": "lead_overlap",
                "eqtl_supported": True,
                "eqtl_lookup_status": "complete",
                "eqtl_lookup_n_hits": 3,
            },
            {
                "trait_id": "trait1",
                "locus_id": "locus2",
                "candidate_rank": 1,
                "gene_id": "ENSG2",
                "gene_name": "GENE2",
                "eqtl_gene_role": "candidate",
                "eqtl_lookup_hit": False,
                "mapping_relation": "nearby_window",
                "eqtl_supported": False,
                "eqtl_lookup_status": "skipped_no_rsid",
                "eqtl_lookup_n_hits": 0,
            },
        ]
    )
    out = build_trait_gene_rollup(gene_candidates)
    assert bool(out.loc[0, "eqtl_lookup_hit"]) is True
    assert out.loc[0, "eqtl_lookup_hit_locus_count"] == 1
    assert out.loc[0, "eqtl_lookup_hit_count"] == 3
    assert out.loc[0, "eqtl_supported_gene_count"] == 1
    assert out.loc[0, "eqtl_supported_locus_count"] == 1
    assert out.loc[0, "eqtl_lookup_complete_locus_count"] == 1
    assert out.loc[0, "eqtl_lookup_skipped_no_rsid_locus_count"] == 1
    assert out.loc[0, "eqtl_total_hit_count"] == 3


def test_build_release_loci_merges_gene_summary_columns():
    build_release_loci = EXPORT_MOD["build_release_loci"]
    loci = pd.DataFrame([{"locus_id": "locus1", "trait_id": "trait1", "lead_rsid": "rs1"}])
    gene_summary = pd.DataFrame(
        [
            {
                "locus_id": "locus1",
                "top_gene_id": "ENSG1",
                "top_gene_name": "GENE1",
                "best_mapping_relation": "lead_overlap",
                "eqtl_lookup_hit": True,
                "eqtl_lookup_hit_gene_count": 1,
                "eqtl_supported": True,
                "eqtl_supported_gene_count": 1,
                "eqtl_study_count": 1,
                "best_eqtl_pvalue": 1e-6,
                "eqtl_lookup_status": "complete",
                "eqtl_lookup_mode": "rsid_api_dataset:v3:QTD000116",
                "eqtl_lookup_n_hits": 2,
            }
        ]
    )
    out = build_release_loci(loci, gene_summary)
    assert out.loc[0, "top_gene_name"] == "GENE1"
    assert bool(out.loc[0, "eqtl_lookup_hit"]) is True
    assert bool(out.loc[0, "eqtl_supported"]) is True
    assert out.loc[0, "eqtl_lookup_status"] == "complete"


def test_build_release_loci_preserves_locus_level_eqtl_summary():
    build_release_loci = EXPORT_MOD["build_release_loci"]
    loci = pd.DataFrame([{"locus_id": "locus1", "trait_id": "trait1", "lead_rsid": "rs1"}])
    gene_summary = pd.DataFrame(
        [
            {
                "locus_id": "locus1",
                "top_gene_id": "ENSG1",
                "top_gene_name": "LOCAL_TOP",
                "best_mapping_relation": "lead_overlap",
                "eqtl_lookup_hit": True,
                "eqtl_lookup_hit_gene_count": 2,
                "eqtl_supported": True,
                "eqtl_supported_gene_count": 1,
                "eqtl_study_count": 2,
                "best_eqtl_pvalue": 1e-7,
                "eqtl_lookup_status": "complete",
                "eqtl_lookup_mode": "rsid_api_dataset:v3:QTD000116",
                "eqtl_lookup_n_hits": 4,
            }
        ]
    )
    out = build_release_loci(loci, gene_summary)
    assert out.loc[0, "top_gene_name"] == "LOCAL_TOP"
    assert bool(out.loc[0, "eqtl_supported"]) is True
    assert out.loc[0, "best_eqtl_pvalue"] == 1e-7
