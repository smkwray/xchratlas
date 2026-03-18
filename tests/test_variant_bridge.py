from __future__ import annotations

import gzip

import pandas as pd

from xatlas.variant_bridge import build_variant_bridge, normalize_variant_key


def test_normalize_variant_key():
    assert normalize_variant_key("chrX", "123", "A", "T") == ("X", 123, "A", "T")


def test_build_variant_bridge_matches_loci(tmp_path):
    manifest = tmp_path / "variant_manifest.txt.bgz"
    rows = [
        "chrom\tpos\tref\talt\trsid\tvarid\n",
        "X\t10\tA\tG\trs10\tX:10_A_G\n",
        "X\t20\tC\tT\tX:20_C_T\tX:20_C_T\n",
        "1\t30\tG\tA\trs30\t1:30_G_A\n",
    ]
    with gzip.open(manifest, "wt", encoding="utf-8") as handle:
        handle.writelines(rows)

    loci = pd.DataFrame(
        {
            "trait_id": ["t1", "t2", "t3"],
            "lead_chr": ["X", "chrX", "X"],
            "lead_pos": [10, 20, 99],
            "lead_ref": ["A", "C", "G"],
            "lead_alt": ["G", "T", "A"],
        }
    )

    annotated, stats = build_variant_bridge(loci, manifest)

    assert stats["unique_keys"] == 3
    assert stats["matched_keys"] == 2
    assert stats["matched_loci"] == 2
    assert stats["matched_rsid_loci"] == 1

    row1 = annotated.loc[annotated["lead_pos"] == 10].iloc[0]
    assert row1["lead_rsid"] == "rs10"
    assert row1["lead_varid"] == "X:10_A_G"

    row2 = annotated.loc[annotated["lead_pos"] == 20].iloc[0]
    assert row2["lead_rsid"] == "X:20_C_T"
    assert row2["lead_varid"] == "X:20_C_T"

    row3 = annotated.loc[annotated["lead_pos"] == 99].iloc[0]
    assert pd.isna(row3["lead_rsid"])
    assert pd.isna(row3["lead_varid"])
