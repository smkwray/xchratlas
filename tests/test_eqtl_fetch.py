import runpy
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
FETCH_MOD = runpy.run_path(str(REPO_ROOT / "scripts/06_fetch_eqtl_region_hits.py"))


class _DummySession:
    pass


def test_lookup_eqtl_by_rsid_aggregates_hits_across_datasets():
    lookup_eqtl_by_rsid = FETCH_MOD["lookup_eqtl_by_rsid"]

    calls = []

    def fake_try_api_request(session, *, url, params, rsid, timeout_seconds):
        calls.append(url)
        if "QTD000001" in url:
            return "hits", pd.DataFrame([{"study_id": "S1", "gene_id": "ENSG1", "pvalue": 1e-8}]), "dataset one"
        if "QTD000002" in url:
            return "hits", pd.DataFrame([{"study_id": "S2", "gene_id": "ENSG2", "pvalue": 5e-6}]), "dataset two"
        raise AssertionError(url)

    lookup_eqtl_by_rsid.__globals__["try_api_request"] = fake_try_api_request

    hits, status, mode, detail = lookup_eqtl_by_rsid(
        _DummySession(),
        rsid="rs123",
        api_base_url="https://example.test",
        timeout_seconds=30.0,
        study_specs=[{"dataset_id": "QTD000001"}, {"dataset_id": "QTD000002"}],
        pvalue_threshold=1e-5,
    )

    assert len(calls) == 2
    assert status == "complete"
    assert mode == "rsid_api_aggregate"
    assert int(hits["eqtl_dataset_id"].nunique()) == 2
    assert set(hits["gene_id"]) == {"ENSG1", "ENSG2"}
    assert "Aggregated" in detail
