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
