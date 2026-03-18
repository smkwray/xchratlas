#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse
import gzip
import io
import shutil
import subprocess
import time
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

from xatlas.constants import DEFAULT_PANUKB_OUTPUT_COLUMNS, X_CHROM_LENGTH_GRCH37
from xatlas.io import ensure_dir, write_tsv
from xatlas.settings import INTERIM_DIR

TABIX_TIMEOUT_SECONDS = 45
REMOTE_TABIX_RETRIES = 3


def public_url_from_manifest_row(row: pd.Series) -> str:
    for field in ["aws_link", "aws_path", "download_url", "url"]:
        if field not in row:
            continue
        value = str(row[field] or "").strip()
        if not value or value.lower() == "nan":
            continue
        if value.startswith("s3://"):
            bucket_and_key = value[len("s3://") :]
            bucket, _, key = bucket_and_key.partition("/")
            if bucket and key:
                return f"https://{bucket}.s3.amazonaws.com/{key}"
        return value
    return ""


def parse_header_and_keep_indices(header_line: str) -> tuple[list[str], list[int]]:
    columns = header_line.rstrip("\n").split("\t")
    keep = [i for i, c in enumerate(columns) if c in DEFAULT_PANUKB_OUTPUT_COLUMNS]
    # Guarantee the key fields even if they were not in the preferred list.
    for required in ["chr", "pos", "ref", "alt"]:
        if required in columns and columns.index(required) not in keep:
            keep.append(columns.index(required))
    keep = sorted(set(keep))
    kept_names = [columns[i] for i in keep]
    return kept_names, keep


def fetch_remote_header_line(url: str) -> str:
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        gz = gzip.GzipFile(fileobj=response.raw)
        text_stream = io.TextIOWrapper(gz, encoding="utf-8")
        header_line = text_stream.readline()
    if not header_line:
        raise RuntimeError("No header received from remote source.")
    return header_line


def reduce_stream_to_chrX(
    line_iter: Iterable[str],
    out_path: Path,
    *,
    already_chrX_filtered: bool = False,
) -> tuple[int, list[str]]:
    row_count = 0
    kept_columns: list[str] = []
    out_path.parent.mkdir(parents=True, exist_ok=True)

    iterator = iter(line_iter)
    header_line = next(iterator, None)
    if header_line is None:
        raise RuntimeError("No header received from source stream.")

    kept_columns, keep_indices = parse_header_and_keep_indices(header_line)

    chr_idx = header_line.rstrip("\n").split("\t").index("chr")
    allowed_chr = {"x", "23", "chrx"}

    with gzip.open(out_path, "wt", encoding="utf-8") as out:
        out.write("\t".join(kept_columns) + "\n")
        for line in iterator:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if not already_chrX_filtered:
                chrom = parts[chr_idx].strip().lower()
                if chrom not in allowed_chr:
                    continue
            kept = [parts[i] if i < len(parts) else "" for i in keep_indices]
            out.write("\t".join(kept) + "\n")
            row_count += 1

    return row_count, kept_columns


def try_remote_tabix(url: str, out_path: Path) -> tuple[str, int, list[str]]:
    if shutil.which("tabix") is None:
        return try_remote_tabix_pysam(url, out_path)

    candidate_regions = [
        "X",
        f"X:1-{X_CHROM_LENGTH_GRCH37}",
        "23",
        f"23:1-{X_CHROM_LENGTH_GRCH37}",
        "chrX",
        f"chrX:1-{X_CHROM_LENGTH_GRCH37}",
    ]

    last_error = None
    for region in candidate_regions:
        cmd = ["tabix", "-h", url, region]
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
                timeout=TABIX_TIMEOUT_SECONDS,
            )
            rows, cols = reduce_stream_to_chrX(
                io.StringIO(result.stdout),
                out_path,
                already_chrX_filtered=True,
            )
            stderr = result.stderr
            rc = result.returncode
            if rc == 0 and rows > 0:
                return f"tabix:{region}", rows, cols
            last_error = RuntimeError(f"tabix returned rc={rc}, rows={rows}, stderr={stderr[:500]}")
        except subprocess.TimeoutExpired:
            if out_path.exists():
                out_path.unlink()
            last_error = RuntimeError(f"tabix timed out after {TABIX_TIMEOUT_SECONDS}s for region {region}")
    raise last_error or RuntimeError("tabix failed for unknown reasons.")


def try_remote_tabix_pysam(url: str, out_path: Path) -> tuple[str, int, list[str]]:
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("Neither system 'tabix' nor Python package 'pysam' is available.") from exc

    candidate_regions = [
        "X",
        f"X:1-{X_CHROM_LENGTH_GRCH37}",
        "23",
        f"23:1-{X_CHROM_LENGTH_GRCH37}",
        "chrX",
        f"chrX:1-{X_CHROM_LENGTH_GRCH37}",
    ]

    tbx = pysam.TabixFile(url)
    try:
        header_line = fetch_remote_header_line(url).rstrip("\n")

        last_error = None
        for region in candidate_regions:
            try:
                iterator = tbx.fetch(region=region)
                line_iter = iter([header_line + "\n", *[line + "\n" for line in iterator]])
                rows, cols = reduce_stream_to_chrX(line_iter, out_path, already_chrX_filtered=True)
                if rows > 0:
                    return f"pysam:{region}", rows, cols
                last_error = RuntimeError(f"pysam returned zero rows for region {region}")
            except Exception as exc:
                if out_path.exists():
                    out_path.unlink()
                last_error = exc
        raise RuntimeError(str(last_error) if last_error else "pysam failed for unknown reasons.")
    finally:
        tbx.close()


def http_stream_chrX(url: str, out_path: Path) -> tuple[str, int, list[str]]:
    with requests.get(url, stream=True, timeout=300) as response:
        response.raise_for_status()
        gz = gzip.GzipFile(fileobj=response.raw)
        text_stream = io.TextIOWrapper(gz, encoding="utf-8")
        rows, cols = reduce_stream_to_chrX(text_stream, out_path, already_chrX_filtered=False)
    return "http_stream_fullfile_filter", rows, cols


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract reduced chrX slices from Pan-UKB per-phenotype files.")
    parser.add_argument("--selected", default=str(INTERIM_DIR / "selected_traits.tsv"))
    parser.add_argument("--outdir", default=str(INTERIM_DIR / "panukb_x"))
    parser.add_argument(
        "--method",
        default="auto",
        choices=["auto", "tabix", "http"],
        help="Preferred extraction method.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    selected = pd.read_csv(args.selected, sep="\t")
    outdir = ensure_dir(args.outdir)
    records = []

    for _, row in selected.iterrows():
        trait_id = row["trait_id"]
        url = public_url_from_manifest_row(row)
        if not url:
            raise ValueError(
                "Selected trait table must contain a usable source path in one of: aws_link, aws_path, download_url, or url."
            )
        out_path = outdir / f"{trait_id}.tsv.gz"

        if out_path.exists() and out_path.stat().st_size == 0:
            out_path.unlink()

        if out_path.exists() and not args.overwrite:
            print(f"[skip] {trait_id} -> {out_path}")
            records.append(
                {
                    "trait_id": trait_id,
                    "source_url": url,
                    "output_path": str(out_path),
                    "method_used": "existing",
                    "rows_written": None,
                    "columns_written": None,
                }
            )
            continue

        print(f"[extract] {trait_id}")
        method_used = None
        rows_written = 0
        columns_written: list[str] = []

        if args.method in {"auto", "tabix"}:
            last_exc = None
            for attempt in range(1, REMOTE_TABIX_RETRIES + 1):
                try:
                    method_used, rows_written, columns_written = try_remote_tabix(url, out_path)
                    last_exc = None
                    break
                except Exception as exc:
                    last_exc = exc
                    if out_path.exists():
                        out_path.unlink()
                    if args.method == "tabix" and attempt == REMOTE_TABIX_RETRIES:
                        raise
                    if attempt < REMOTE_TABIX_RETRIES:
                        print(f"[warn] remote tabix attempt {attempt}/{REMOTE_TABIX_RETRIES} failed for {trait_id}: {exc}")
                        time.sleep(2.0 * attempt)
            if method_used is None and last_exc is not None:
                print(f"[warn] tabix failed for {trait_id}: {last_exc}")

        if method_used is None:
            method_used, rows_written, columns_written = http_stream_chrX(url, out_path)

        records.append(
            {
                "trait_id": trait_id,
                "source_url": url,
                "output_path": str(out_path),
                "method_used": method_used,
                "rows_written": rows_written,
                "columns_written": ",".join(columns_written),
            }
        )
        print(f"[done] {trait_id}: {rows_written:,} rows via {method_used}")

    manifest = pd.DataFrame(records)
    write_tsv(manifest, INTERIM_DIR / "x_extract_manifest.tsv")
    print(f"[done] chrX extraction manifest written to {INTERIM_DIR / 'x_extract_manifest.tsv'}")


if __name__ == "__main__":
    main()
