import pandas as pd
import pytest

from xatlas.scoring import lookup_counts

from xatlas.scoring import compute_trait_score, coverage_grade, effective_hq_case_count, summarize_eqtl_lookup_note


def test_summarize_eqtl_lookup_note_handles_partial_rsid_coverage():
    rows = pd.DataFrame(
        [
            {
                "locus_id": "locus1",
                "eqtl_lookup_status": "complete",
                "eqtl_lookup_mode": "rsid_api_dataset:v3:QTD000116",
                "eqtl_lookup_n_hits": 3,
            },
            {
                "locus_id": "locus2",
                "eqtl_lookup_status": "skipped_no_rsid",
                "eqtl_lookup_mode": "rsid_api",
                "eqtl_lookup_n_hits": 0,
            },
        ]
    )
    assert (
        summarize_eqtl_lookup_note(rows)
        == "eQTL follow-up completed for rsID-addressable loci; some loci lacked a usable rsID."
    )


def test_summarize_eqtl_lookup_note_keeps_complete_no_support_message():
    rows = pd.DataFrame(
        [
            {
                "locus_id": "locus1",
                "eqtl_lookup_status": "no_hits",
                "eqtl_lookup_mode": "rsid_api_dataset:v3:QTD000116",
                "eqtl_lookup_n_hits": 0,
            },
            {
                "locus_id": "locus2",
                "eqtl_lookup_status": "complete",
                "eqtl_lookup_mode": "rsid_api_dataset:v3:QTD000121",
                "eqtl_lookup_n_hits": 0,
            },
        ]
    )
    assert summarize_eqtl_lookup_note(rows) == "eQTL follow-up: no support observed in lookup summary."


def test_summarize_eqtl_lookup_note_reports_lookup_hit_only_state():
    rows = pd.DataFrame(
        [
            {
                "locus_id": "locus1",
                "eqtl_lookup_status": "complete",
                "eqtl_lookup_mode": "rsid_api_aggregate",
                "eqtl_lookup_n_hits": 5,
                "eqtl_supported": False,
            }
        ]
    )
    assert (
        summarize_eqtl_lookup_note(rows)
        == "eQTL follow-up returned lookup-hit associations, but none met the strict candidate-gene support rule."
    )


def test_effective_hq_case_count_uses_requested_stratum_sample_size():
    row = {
        "pheno_sex": "females",
        "n_cases_hq_cohort_both_sexes": 100000,
        "n_cases_hq_cohort_females": 12000,
        "n_cases_hq_cohort_males": 8000,
    }
    assert effective_hq_case_count(row) == 12000


def test_coverage_grade_uses_stratum_specific_case_count():
    row = {
        "pheno_sex": "males",
        "in_max_independent_set": True,
        "num_pops_pass_qc": 1,
        "n_cases_hq_cohort_both_sexes": 75000,
        "n_cases_hq_cohort_males": 5000,
    }
    assert coverage_grade(row) == "C"


def test_compute_trait_score_adds_sex_specific_confidence_notes():
    score = compute_trait_score(
        {
            "pheno_sex": "males",
            "in_max_independent_set": True,
            "num_pops_pass_qc": 2,
            "n_cases_hq_cohort_males": 60000,
        },
        {
            "n_loci": 1,
            "lead_neglog10_pvalue": 9.5,
            "any_nonpar_locus": True,
            "par_only_signal": False,
            "pvalue_scale": "neglog10",
        },
        gene_row={"best_mapping_relation": "lead_overlap", "eqtl_supported": False},
    )
    assert "male-only" in score["confidence_notes"]

@pytest.mark.parametrize(("observations", "expected"), [
    ([], (0, 0, 0, 0)),
    ([("complete", 0)], (1, 1, 0, 0)),
    ([("skipped_no_rsid", 0)], (0, 0, 1, 0)),
    ([("complete", 3), ("skipped_no_rsid", 0)], (1, 0, 1, 1)),
])
def test_lookup_assessment_counts(observations, expected):
    counts = lookup_counts([{"eqtl_lookup_status": status, "eqtl_lookup_n_hits": hits}
                            for status, hits in observations])
    assert tuple(counts[key] for key in ("eqtl_lookup_complete_locus_count",
        "eqtl_lookup_no_hit_locus_count", "eqtl_lookup_skipped_no_rsid_locus_count",
        "eqtl_lookup_hit_locus_count")) == expected
