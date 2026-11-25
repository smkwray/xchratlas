import math
import runpy
import sys
from pathlib import Path

import pytest
import pandas as pd

from xatlas.loci import call_loci_from_dataframe, infer_pvalue_scale, to_raw_pvalue


def test_infer_scale_neglog10():
    s = pd.Series([8.0, 10.2, 7.5])
    assert infer_pvalue_scale(s, "neglog10_pval_meta") == "neglog10"


def test_infer_scale_raw():
    s = pd.Series([0.5, 1e-5, 0.01])
    assert infer_pvalue_scale(s, "pval_meta") == "raw"


def test_to_raw_pvalue_ln():
    s = pd.Series([-1.0, -10.0])
    p, scale = to_raw_pvalue(s, "log_pval_meta")
    assert scale == "ln"
    assert (p > 0).all()
    assert (p <= 1).all()


def test_call_loci_simple():
    df = pd.DataFrame(
        {
            "chr": ["X", "X", "X", "X"],
            "pos": [1_000_000, 1_200_000, 10_000_000, 10_100_000],
            "ref": ["A", "A", "G", "G"],
            "alt": ["C", "T", "A", "T"],
            "neglog10_pval_meta_hq": [9.0, 8.0, 12.0, 2.0],
            "beta_meta_hq": [0.1, 0.2, -0.3, 0.4],
            "se_meta_hq": [0.01, 0.02, 0.03, 0.04],
        }
    )
    loci, summary = call_loci_from_dataframe(df, trait_id="test_trait", window_bp=500_000)
    assert summary["n_loci"] == 2
    assert len(loci) == 2


@pytest.mark.parametrize(("column", "values", "expected"), [
    ("neglog10_pval_meta", [310.0, 400.0], 400.0),
    ("neglog10_pval_meta", [10, 20], 20.0),
    ("log_pval_meta", [-310 * math.log(10), -400 * math.log(10)], 400.0),
    ("pval_meta", [1e-305, 1e-310], 310.0),
])
def test_log_space_lead(column, values, expected):
    frame = pd.DataFrame({"chr": ["X", "X"], "pos": [100, 200], column: values})
    for rows in (frame, frame.iloc[::-1]):
        loci, summary = call_loci_from_dataframe(rows, trait_id="test")
        assert loci.iloc[0]["lead_pos"] == 200
        assert summary["lead_neglog10_pvalue"] == pytest.approx(expected)


@pytest.mark.parametrize("extract", ["missing", "invalid", "nonsignificant"])
def test_extract_result(tmp_path, monkeypatch, extract):
    script = Path(__file__).resolve().parents[1] / "scripts/05_call_x_loci.py"
    selected = tmp_path / "selected.tsv"
    selected.write_text("trait_id\ntrait1\n")
    summary = tmp_path / "summary.tsv"
    monkeypatch.setattr(sys, "argv", [str(script), "--extract-dir", str(tmp_path),
        "--selected", str(selected), "--out-loci", str(tmp_path / "loci.tsv.gz"),
        "--out-trait-summary", str(summary)])
    frame = pd.DataFrame({"chr": ["X"], "pos": [123], "pval_meta": [0.5]})
    if extract == "invalid":
        frame = frame.drop(columns="pval_meta")
    if extract != "missing":
        frame.to_csv(tmp_path / "trait1.tsv.gz", sep="\t", index=False)
    if extract == "nonsignificant":
        runpy.run_path(str(script), run_name="__main__")
        assert pd.read_csv(summary, sep="\t").iloc[0]["n_loci"] == 0
    else:
        with pytest.raises((FileNotFoundError, ValueError), match="trait1"):
            runpy.run_path(str(script), run_name="__main__")
        assert not summary.exists()
