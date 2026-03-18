#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse
import io
import re
from pathlib import Path

import pandas as pd
import requests

from xatlas.constants import (
    EQTL_CHRX_GENOTYPES_URL,
    EQTL_TABIX_IMPORTED_PATHS_URL,
    EQTL_TABIX_PATHS_URL,
)
from xatlas.io import choose_column, ensure_dir, write_tsv
from xatlas.settings import INTERIM_DIR, RAW_DIR


def download_text(url: str) -> str:
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    return response.text


def robust_read_tsv(url: str) -> pd.DataFrame:
    text = download_text(url)
    df = pd.read_csv(io.StringIO(text), sep="\t")
    if df.shape[1] > 1:
        return df
    return pd.read_csv(io.StringIO(text), sep=r"\s+", engine="python")


def read_chrX_support(url: str) -> pd.DataFrame:
    text = download_text(url)
    try:
        df = pd.read_csv(io.StringIO(text), sep="\t")
        if df.shape[1] >= 2 and "study_id" in [c.lower() for c in df.columns]:
            return df
    except Exception:
        pass

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        return pd.DataFrame(columns=["study_id", "chrx_status", "comment"])

    rows = []
    for line in lines[1:]:
        parts = line.split(maxsplit=2)
        rows.append(
            {
                "study_id": parts[0] if len(parts) > 0 else "",
                "chrx_status": parts[1] if len(parts) > 1 else "",
                "comment": parts[2] if len(parts) > 2 else "",
            }
        )
    return pd.DataFrame(rows)


def annotate_paths(df: pd.DataFrame, source_kind: str) -> pd.DataFrame:
    out = df.copy()
    out["source_kind"] = source_kind
    return out


def canonicalize_study_id(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    aliases = {
        "GTEx_V8": "GTEx",
    }
    if text in aliases:
        return aliases[text]
    text = re.sub(r"_V\d+$", "", text)
    return text


def first_non_empty_value(*values: object) -> str:
    for value in values:
        text = str(value or "").strip()
        if text and text.lower() != "nan":
            return text
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare a small eQTL Catalogue chrX index.")
    parser.add_argument("--outdir", default=str(RAW_DIR / "eqtl_catalogue"))
    parser.add_argument("--priority-studies", default="config/eqtl_priority_studies.csv")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)

    print("[download] eQTL chrX support metadata")
    chrx = read_chrX_support(EQTL_CHRX_GENOTYPES_URL)
    if "study_id" not in chrx.columns:
        # normalise first column name if needed
        chrx.columns = [str(c).strip().lower().replace(" ", "_") for c in chrx.columns]
        if "study_id" not in chrx.columns and len(chrx.columns) >= 1:
            chrx = chrx.rename(columns={chrx.columns[0]: "study_id"})
    if "chrx_status" not in chrx.columns:
        # heuristically rename the first non-study column
        for col in chrx.columns:
            if col != "study_id":
                chrx = chrx.rename(columns={col: "chrx_status"})
                break
    write_tsv(chrx, outdir / "chrX_genotypes.tsv")

    print("[download] eQTL tabix path metadata")
    uniform = annotate_paths(robust_read_tsv(EQTL_TABIX_PATHS_URL), "uniform")
    imported = annotate_paths(robust_read_tsv(EQTL_TABIX_IMPORTED_PATHS_URL), "imported")
    write_tsv(uniform, outdir / "tabix_ftp_paths.tsv")
    write_tsv(imported, outdir / "tabix_ftp_paths_imported.tsv")

    paths = pd.concat([uniform, imported], ignore_index=True, sort=False)

    study_accession_col = choose_column(list(paths.columns), exact=["study_id"], contains_all=["study", "id"])
    study_label_col = choose_column(list(paths.columns), exact=["study_label", "study"], contains_any=["study_label", "study"])
    dataset_col = choose_column(list(paths.columns), exact=["dataset_id"], contains_all=["dataset", "id"])
    ftp_col = choose_column(list(paths.columns), exact=["ftp_path"], contains_all=["ftp", "path"]) or choose_column(
        list(paths.columns),
        contains_any=["ftp", "url", "path"],
    )
    qtl_group_col = choose_column(list(paths.columns), exact=["qtl_group", "sample_group"], contains_any=["qtl_group", "sample_group"])
    quant_method_col = choose_column(list(paths.columns), exact=["quant_method"], contains_all=["quant", "method"])
    tissue_col = choose_column(list(paths.columns), contains_any=["tissue", "cell_type", "celltype"])

    if study_accession_col is None and study_label_col is None:
        raise RuntimeError("Could not find a study column in eQTL path metadata.")
    if ftp_col is None:
        raise RuntimeError("Could not find an FTP/path column in eQTL path metadata.")

    chrx_simple = chrx.copy()
    if "study_id" not in chrx_simple.columns:
        raise RuntimeError("chrX support table is missing a study_id column.")
    if "chrx_status" not in chrx_simple.columns:
        raise RuntimeError("chrX support table is missing a chrX status column.")
    chrx_simple["study_id_canonical"] = chrx_simple["study_id"].map(canonicalize_study_id)
    paths["study_accession"] = (
        paths[study_accession_col].astype("string") if study_accession_col else pd.Series(pd.NA, index=paths.index, dtype="string")
    )
    paths["study_label_raw"] = paths.apply(
        lambda row: first_non_empty_value(
            row.get(study_label_col) if study_label_col else "",
            row.get(study_accession_col) if study_accession_col else "",
        ),
        axis=1,
    )
    paths["study_id_canonical"] = paths["study_label_raw"].map(canonicalize_study_id)
    paths["qtl_group"] = paths.apply(
        lambda row: first_non_empty_value(
            row.get(qtl_group_col) if qtl_group_col else "",
            row.get(dataset_col) if dataset_col else "",
        ),
        axis=1,
    )

    merged = paths.merge(chrx_simple, on="study_id_canonical", how="left", suffixes=("", "_support"))
    merged["study_id"] = merged["study_id_canonical"]
    merged["chrX_supported"] = merged["chrx_status"].fillna("").astype(str).str.lower().isin({"yes", "nonpar"})
    merged["ftp_path_inferred"] = merged[ftp_col]

    if Path(args.priority_studies).exists():
        priority = pd.read_csv(args.priority_studies)
        merged = merged.merge(priority, on="study_id", how="left")
    else:
        merged["priority"] = pd.NA

    keep_cols = [
        "study_id",
        "study_accession",
        "source_kind",
        "chrx_status",
        "chrX_supported",
        "ftp_path_inferred",
        "priority",
        "qtl_group",
    ]
    for col in [study_label_col, dataset_col, quant_method_col, tissue_col]:
        if col and col not in keep_cols:
            keep_cols.append(col)
    keep_cols += [c for c in ["comment"] if c in merged.columns]
    keep_cols = [c for c in dict.fromkeys(keep_cols) if c in merged.columns]

    eqtl_index = merged[keep_cols].drop_duplicates().sort_values(
        ["priority", "study_id"], na_position="last"
    )
    write_tsv(eqtl_index, INTERIM_DIR / "eqtl_chrX_index.tsv")

    print(f"[done] wrote {len(eqtl_index):,} eQTL rows to {INTERIM_DIR / 'eqtl_chrX_index.tsv'}")


if __name__ == "__main__":
    main()
