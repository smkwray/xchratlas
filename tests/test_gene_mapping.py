import runpy
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
MAP_GENES_MOD = runpy.run_path(str(REPO_ROOT / "scripts/07_map_loci_to_genes.py"))


def test_append_missing_eqtl_hit_genes_appends_external_eqtl_gene_after_structural_candidates():
    append_missing_eqtl_hit_genes = MAP_GENES_MOD["append_missing_eqtl_hit_genes"]

    genes = pd.DataFrame(
        [
            {
                "gene_id": "ENSG000001.1",
                "gene_id_clean": "ENSG000001",
                "gene_name": "LOCAL1",
                "gene_biotype": "protein_coding",
                "start": 100,
                "end": 200,
            },
            {
                "gene_id": "ENSG000002.1",
                "gene_id_clean": "ENSG000002",
                "gene_name": "LOCAL2",
                "gene_biotype": "protein_coding",
                "start": 250,
                "end": 350,
            },
            {
                "gene_id": "ENSG000003.1",
                "gene_id_clean": "ENSG000003",
                "gene_name": "EQTL1",
                "gene_biotype": "protein_coding",
                "start": 1200,
                "end": 1300,
            },
        ]
    )
    base_subset = genes.iloc[:2].copy()
    base_subset["distance_to_lead_bp"] = [0, 50]
    base_subset["mapping_relation"] = ["lead_overlap", "locus_overlap"]
    eqtl_support = {"ENSG000003": {"eqtl_lookup_hit": True, "eqtl_supported": True, "best_eqtl_pvalue": 1e-6}}

    out = append_missing_eqtl_hit_genes(
        base_subset,
        genes=genes,
        eqtl_support=eqtl_support,
        lead_pos=150,
        locus_start=100,
        locus_end=350,
        flank_bp=100,
    )

    assert list(out["gene_id_clean"]) == ["ENSG000001", "ENSG000002", "ENSG000003"]
    assert out.iloc[2]["mapping_relation"] == "eqtl_hit_gene"
    assert out.iloc[2]["eqtl_gene_role"] == "followup"


def test_append_missing_eqtl_hit_genes_ignores_eqtl_genes_already_present():
    append_missing_eqtl_hit_genes = MAP_GENES_MOD["append_missing_eqtl_hit_genes"]

    genes = pd.DataFrame(
        [
            {
                "gene_id": "ENSG000001.1",
                "gene_id_clean": "ENSG000001",
                "gene_name": "LOCAL1",
                "gene_biotype": "protein_coding",
                "start": 100,
                "end": 200,
            }
        ]
    )
    base_subset = genes.copy()
    base_subset["distance_to_lead_bp"] = [0]
    base_subset["mapping_relation"] = ["lead_overlap"]
    eqtl_support = {"ENSG000001": {"eqtl_lookup_hit": True, "eqtl_supported": True, "best_eqtl_pvalue": 1e-6}}

    out = append_missing_eqtl_hit_genes(
        base_subset,
        genes=genes,
        eqtl_support=eqtl_support,
        lead_pos=150,
        locus_start=100,
        locus_end=200,
        flank_bp=100,
    )

    assert len(out) == 1


def test_append_missing_eqtl_hit_genes_keeps_lookup_only_external_genes_as_followup():
    append_missing_eqtl_hit_genes = MAP_GENES_MOD["append_missing_eqtl_hit_genes"]

    genes = pd.DataFrame(
        [
            {
                "gene_id": "ENSG000001.1",
                "gene_id_clean": "ENSG000001",
                "gene_name": "LOCAL1",
                "gene_biotype": "protein_coding",
                "start": 100,
                "end": 200,
            },
            {
                "gene_id": "ENSG000003.1",
                "gene_id_clean": "ENSG000003",
                "gene_name": "WEAK_EQTL",
                "gene_biotype": "protein_coding",
                "start": 1200,
                "end": 1300,
            },
        ]
    )
    base_subset = genes.iloc[:1].copy()
    base_subset["distance_to_lead_bp"] = [0]
    base_subset["mapping_relation"] = ["lead_overlap"]
    eqtl_support = {"ENSG000003": {"eqtl_lookup_hit": True, "eqtl_supported": False, "best_eqtl_pvalue": 0.4}}

    out = append_missing_eqtl_hit_genes(
        base_subset,
        genes=genes,
        eqtl_support=eqtl_support,
        lead_pos=150,
        locus_start=100,
        locus_end=200,
        flank_bp=100,
    )

    assert list(out["gene_id_clean"]) == ["ENSG000001", "ENSG000003"]
    assert out.iloc[1]["eqtl_gene_role"] == "followup"


def test_gene_distance_to_lead_recomputes_cleanly_for_each_gene():
    gene_distance_to_lead = MAP_GENES_MOD["gene_distance_to_lead"]

    gene = pd.Series({"start": 51804923, "end": 51812368})
    assert gene_distance_to_lead(gene, 51705012) == 99911
