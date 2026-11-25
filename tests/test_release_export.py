import runpy
import sys

import pytest
from xatlas.scoring import compute_trait_score
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


def test_returned_records_survive_an_unmapped_lookup():
    loci = pd.DataFrame([{"locus_id": "locus1", "trait_id": "trait1"}])
    summary = pd.DataFrame([{"locus_id": "locus1", "eqtl_lookup_hit": False,
                             "eqtl_lookup_n_hits": 3, "eqtl_lookup_hit_gene_count": 0}])
    released = EXPORT_MOD["build_release_loci"](loci, summary)
    assert bool(released.loc[0, "eqtl_lookup_hit"])
    assert released.loc[0, "eqtl_lookup_hit_gene_count"] == 0


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


def test_release_export_uses_source_hq_count(tmp_path: Path):
    selected = pd.DataFrame(
        [
            {
                "trait_id": "trait1",
                "query_id": "trait_one",
                "description": "Trait One",
                "domain": "alpha",
                "pheno_sex": "females",
                "num_pops_pass_qc": 2,
                "in_max_independent_set": True,
                "n_cases_hq_cohort_females": 12000,
            }
        ]
    )
    selected_path = tmp_path / "selected.tsv"
    selected.to_csv(selected_path, sep="\t", index=False)

    trait_summary = pd.DataFrame(
        [
            {
                "trait_id": "trait1",
                "n_loci": 0,
                "lead_neglog10_pvalue": None,
                "pvalue_scale": "neglog10",
                "par_only_signal": False,
                "any_nonpar_locus": False,
            }
        ]
    )
    trait_summary_path = tmp_path / "trait_summary.tsv"
    trait_summary.to_csv(trait_summary_path, sep="\t", index=False)

    outdir = tmp_path / "release"
    export_main = EXPORT_MOD["main"]

    import sys

    argv = sys.argv[:]
    try:
        sys.argv = [
            "08_export_release_tables.py",
            "--selected",
            str(selected_path),
            "--trait-summary",
            str(trait_summary_path),
            "--outdir",
            str(outdir),
            "--loci",
            str(tmp_path / "missing_loci.tsv.gz"),
            "--gene-summary",
            str(tmp_path / "missing_gene_summary.tsv"),
            "--gene-candidates",
            str(tmp_path / "missing_gene_candidates.tsv.gz"),
        ]
        export_main()
    finally:
        sys.argv = argv

    trait_scores = pd.read_csv(outdir / "trait_scores.tsv", sep="\t")
    assert trait_scores.loc[0, "effective_hq_case_count"] == 12000

def test_mapping_bonus_is_independent_of_locus_order():
    rows = pd.DataFrame([
        dict(trait_id="t", locus_id=locus, gene_id=gene, gene_name=gene, candidate_rank=1,
             mapping_relation=relation, eqtl_gene_role=role, eqtl_lookup_status="complete",
             eqtl_lookup_n_hits=3, eqtl_lookup_hit=False)
        for locus, gene, relation, role in [
            ("a", "g1", "nearby_window", "candidate"),
            ("b", "g2", "locus_overlap", "candidate"),
            ("c", "g3", "lead_overlap", "followup"),
        ]
    ])
    for overlap in ("locus_overlap", "lead_overlap"):
        rows.loc[1, "mapping_relation"] = overlap
        rollups = [EXPORT_MOD["build_trait_gene_rollup"](order).iloc[0].to_dict()
                   for order in (rows, rows.iloc[::-1])]
        assert rollups[0] == rollups[1]
        assert rollups[0]["best_mapping_relation"] == overlap
        assert rollups[0]["eqtl_lookup_hit_locus_count"] == 3
        assert compute_trait_score({}, {"n_loci": 3}, gene_row=rollups[0])["x_evidence_score"] == 32


@pytest.mark.parametrize("summary", ["absent", "unsuccessful"])
def test_missing_successful_summary(tmp_path, monkeypatch, summary):
    selected = tmp_path / "selected.tsv"
    selected.write_text("trait_id\ntest_trait\n")
    summary_path = tmp_path / "summary.tsv"
    if summary == "unsuccessful":
        summary_path.write_text("trait_id\tn_loci\tprocessing_note\ntest_trait\t0\tmissing extract\n")
    monkeypatch.setattr(sys, "argv", ["export", "--selected", str(selected),
        "--trait-summary", str(summary_path), "--outdir", str(tmp_path / "release")])
    with pytest.raises(ValueError, match="test_trait"):
        EXPORT_MOD["main"]()


def test_no_phenotype_match_is_ineligible():
    selector = runpy.run_path(str(REPO_ROOT / "scripts/01_select_seed_traits.py"))
    row = pd.Series(dict(description="Age at menarche", _search_blob="age at menarche",
                         pheno_sex="females", num_pops_pass_qc=4, in_max_independent_set=True))
    assert selector["query_score"](row, "Standing height", "continuous", "females") == 0
