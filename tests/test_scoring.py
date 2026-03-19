import pandas as pd

from xatlas.scoring import summarize_eqtl_lookup_note


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
