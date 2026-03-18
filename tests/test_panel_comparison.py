import runpy
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
COMPARE_MOD = runpy.run_path(str(REPO_ROOT / "scripts/09_compare_release_panels.py"))


def test_build_panel_release_comparison_computes_presence_and_deltas():
    build_panel_release_comparison = COMPARE_MOD["build_panel_release_comparison"]

    small = pd.DataFrame(
        [
            {
                "trait_id": "trait_a",
                "query_id": "a",
                "domain": "blood",
                "description": "Trait A",
                "n_loci": 1,
                "x_evidence_score": 2.5,
                "eqtl_supported": True,
                "eqtl_supported_locus_count": 1,
                "eqtl_total_hit_count": 4,
                "coverage_grade": "A",
            }
        ]
    )
    expanded = pd.DataFrame(
        [
            {
                "trait_id": "trait_a",
                "query_id": "a",
                "domain": "blood",
                "description": "Trait A",
                "n_loci": 3,
                "x_evidence_score": 4.0,
                "eqtl_supported": True,
                "eqtl_supported_locus_count": 2,
                "eqtl_total_hit_count": 9,
                "coverage_grade": "A",
            },
            {
                "trait_id": "trait_b",
                "query_id": "b",
                "domain": "immune",
                "description": "Trait B",
                "n_loci": 2,
                "x_evidence_score": 3.0,
                "eqtl_supported": "false",
                "eqtl_supported_locus_count": 0,
                "eqtl_total_hit_count": 0,
                "coverage_grade": "B",
            },
        ]
    )

    out = build_panel_release_comparison(small, expanded)

    trait_a = out.loc[out["trait_id"] == "trait_a"].iloc[0]
    assert bool(trait_a["in_small_release"]) is True
    assert bool(trait_a["in_panel_expanded"]) is True
    assert trait_a["delta_n_loci"] == 2
    assert trait_a["delta_x_evidence_score"] == 1.5
    assert trait_a["delta_eqtl_supported_locus_count"] == 1
    assert trait_a["delta_eqtl_total_hit_count"] == 5
    assert trait_a["eqtl_supported_delta"] == 0

    trait_b = out.loc[out["trait_id"] == "trait_b"].iloc[0]
    assert bool(trait_b["in_small_release"]) is False
    assert bool(trait_b["in_panel_expanded"]) is True
    assert trait_b["n_loci_small"] == 0
    assert trait_b["n_loci_expanded"] == 2
    assert trait_b["eqtl_supported_delta"] == 0


def test_build_panel_release_comparison_normalizes_bool_columns():
    build_panel_release_comparison = COMPARE_MOD["build_panel_release_comparison"]

    small = pd.DataFrame(
        [
            {
                "trait_id": "trait_a",
                "query_id": "a",
                "domain": "blood",
                "description": "Trait A",
                "n_loci": 1,
                "x_evidence_score": 2.0,
                "eqtl_supported": "yes",
                "eqtl_supported_locus_count": 1,
                "eqtl_total_hit_count": 3,
                "coverage_grade": "A",
            }
        ]
    )
    expanded = pd.DataFrame(
        [
            {
                "trait_id": "trait_a",
                "query_id": "a",
                "domain": "blood",
                "description": "Trait A",
                "n_loci": 1,
                "x_evidence_score": 2.0,
                "eqtl_supported": "",
                "eqtl_supported_locus_count": 0,
                "eqtl_total_hit_count": 0,
                "coverage_grade": "A",
            }
        ]
    )

    out = build_panel_release_comparison(small, expanded)
    trait_a = out.loc[out["trait_id"] == "trait_a"].iloc[0]
    assert bool(trait_a["eqtl_supported_small"]) is True
    assert bool(trait_a["eqtl_supported_expanded"]) is False
    assert trait_a["eqtl_supported_delta"] == -1
