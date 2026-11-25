#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse
import time
from pathlib import Path

import pandas as pd
import requests

from xatlas.constants import ENSEMBL_GRCH37_REST, X_CHROM_LENGTH_GRCH37
from xatlas.io import write_tsv
from xatlas.region import classify_interval_region
from xatlas.settings import INTERIM_DIR

CHUNK_SIZE = 5_000_000


def fetch_window(session: requests.Session, start: int, end: int) -> list[dict]:
    url = f"{ENSEMBL_GRCH37_REST}/overlap/region/human/X:{start}-{end}"
    response = session.get(
        url,
        params={"feature": "gene"},
        headers={"Content-Type": "application/json"},
        timeout=120,
    )
    response.raise_for_status()
    return response.json()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a GRCh37 chromosome X gene catalog from Ensembl REST.")
    parser.add_argument("--out", default=str(INTERIM_DIR / "ensembl_x_genes.tsv"))
    parser.add_argument("--sleep-seconds", type=float, default=0.1)
    args = parser.parse_args()

    records = []
    session = requests.Session()

    for start in range(1, X_CHROM_LENGTH_GRCH37 + 1, CHUNK_SIZE):
        end = min(X_CHROM_LENGTH_GRCH37, start + CHUNK_SIZE - 1)
        print(f"[fetch] X:{start}-{end}")
        features = fetch_window(session, start, end)
        for feat in features:
            gene_id = feat.get("id")
            gene_start = int(feat["start"])
            gene_end = int(feat["end"])
            records.append(
                {
                    "gene_id": gene_id,
                    "gene_name": feat.get("external_name"),
                    "gene_biotype": feat.get("biotype"),
                    "chromosome": feat.get("seq_region_name"),
                    "start": gene_start,
                    "end": gene_end,
                    "strand": feat.get("strand"),
                    "description": feat.get("description"),
                    "x_region": classify_interval_region(gene_start, gene_end),
                }
            )
        time.sleep(args.sleep_seconds)

    genes = pd.DataFrame(records).drop_duplicates(subset=["gene_id"]).sort_values(["start", "end"])
    write_tsv(genes, Path(args.out))
    print(f"[done] wrote {len(genes):,} genes to {args.out}")


if __name__ == "__main__":
    main()
