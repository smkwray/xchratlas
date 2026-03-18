from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Any

import pandas as pd

from xatlas.eqtl import normalize_rsid


def normalize_variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> tuple[str, int, str, str]:
    chrom_text = str(chrom).replace("chr", "").strip().upper()
    return chrom_text, int(pos), str(ref), str(alt)


def build_variant_bridge(
    loci_df: pd.DataFrame,
    variant_manifest_path: str | Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    path = Path(variant_manifest_path)
    if loci_df.empty:
        return loci_df.copy(), {
            "variant_manifest_path": str(path),
            "unique_keys": 0,
            "matched_keys": 0,
            "matched_loci": 0,
        }
    if not path.exists():
        raise FileNotFoundError(f"Variant manifest not found: {path}")

    required = ["lead_chr", "lead_pos", "lead_ref", "lead_alt"]
    missing = [col for col in required if col not in loci_df.columns]
    if missing:
        raise ValueError(f"Locus table is missing required columns: {', '.join(missing)}")

    loci = loci_df.copy()
    loci["_variant_key"] = [
        normalize_variant_key(chrom, pos, ref, alt)
        for chrom, pos, ref, alt in loci[required].itertuples(index=False, name=None)
    ]
    wanted_keys = set(loci["_variant_key"])
    wanted_chroms = {chrom for chrom, _, _, _ in wanted_keys}

    found: dict[tuple[str, int, str, str], dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = str(row.get("chrom", "")).replace("chr", "").strip().upper()
            if chrom not in wanted_chroms:
                continue
            key = normalize_variant_key(chrom, row["pos"], row["ref"], row["alt"])
            if key not in wanted_keys or key in found:
                continue
            found[key] = {
                "lead_rsid": str(row.get("rsid", "") or ""),
                "lead_varid": str(row.get("varid", "") or ""),
            }
            if len(found) == len(wanted_keys):
                break

    annotations = pd.DataFrame(
        [{"_variant_key": key, **value} for key, value in found.items()]
    )
    if annotations.empty:
        loci["lead_rsid"] = loci.get("lead_rsid", pd.Series(pd.NA, index=loci.index))
        loci["lead_varid"] = loci.get("lead_varid", pd.Series(pd.NA, index=loci.index))
    else:
        loci = loci.drop(columns=[c for c in ["lead_rsid", "lead_varid"] if c in loci.columns]).merge(
            annotations,
            on="_variant_key",
            how="left",
        )

    loci = loci.drop(columns=["_variant_key"])
    loci["lead_rsid"] = loci["lead_rsid"].replace({"": pd.NA})
    loci["lead_varid"] = loci["lead_varid"].replace({"": pd.NA})

    stats = {
        "variant_manifest_path": str(path),
        "unique_keys": len(wanted_keys),
        "matched_keys": len(found),
        "matched_loci": int(loci["lead_varid"].notna().sum()),
        "matched_rsid_loci": int(loci["lead_rsid"].map(normalize_rsid).notna().sum()),
    }
    return loci, stats
