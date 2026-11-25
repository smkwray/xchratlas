#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse
import gzip
import shutil
from pathlib import Path

import requests

from xatlas.constants import (
    PANUKB_PHENOTYPE_MANIFEST_URL,
    PANUKB_VARIANT_MANIFEST_URL,
)
from xatlas.io import ensure_dir
from xatlas.settings import RAW_DIR


def download_file(url: str, destination: Path, overwrite: bool = False) -> Path:
    if destination.exists() and not overwrite:
        print(f"[skip] {destination} already exists")
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f"[download] {url} -> {destination}")
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with open(destination, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    return destination


def maybe_decompress_bgz(path: Path) -> Path:
    if path.suffix not in {".bgz", ".gz"}:
        return path
    plain_name = path.name
    for suffix in [".bgz", ".gz"]:
        if plain_name.endswith(suffix):
            plain_name = plain_name[: -len(suffix)]
            break
    out = path.with_name(plain_name)
    if out.exists():
        print(f"[skip] {out} already exists")
        return out
    print(f"[decompress] {path} -> {out}")
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Download Pan-UKB manifests.")
    parser.add_argument("--outdir", default=str(RAW_DIR / "panukb"), help="Output directory.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    parser.add_argument(
        "--with-variant-manifest",
        action="store_true",
        help="Also download the larger Pan-UKB variant manifest.",
    )
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    phenotype_bgz = download_file(
        PANUKB_PHENOTYPE_MANIFEST_URL,
        outdir / "phenotype_manifest.tsv.bgz",
        overwrite=args.overwrite,
    )
    maybe_decompress_bgz(phenotype_bgz)

    if args.with_variant_manifest:
        download_file(
            PANUKB_VARIANT_MANIFEST_URL,
            outdir / "full_variant_qc_metrics.txt.bgz",
            overwrite=args.overwrite,
        )

    print("[done] Pan-UKB metadata fetch complete.")


if __name__ == "__main__":
    main()
