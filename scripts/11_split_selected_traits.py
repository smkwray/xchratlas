#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import argparse

import pandas as pd

from xatlas.io import ensure_dir, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser(description="Split a selected-traits table into balanced shards.")
    parser.add_argument("--selected", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--num-shards", type=int, default=4)
    parser.add_argument("--weight-column", default="size_in_bytes")
    args = parser.parse_args()

    selected = pd.read_csv(args.selected, sep="\t")
    outdir = ensure_dir(args.outdir)

    weights = pd.to_numeric(selected.get(args.weight_column), errors="coerce").fillna(1.0)
    selected = selected.assign(_weight=weights).sort_values(["_weight", "selection_rank"], ascending=[False, True]).reset_index(drop=True)

    shards: list[list[dict]] = [[] for _ in range(args.num_shards)]
    shard_weights = [0.0] * args.num_shards

    for _, row in selected.iterrows():
        shard_idx = min(range(args.num_shards), key=lambda i: shard_weights[i])
        shards[shard_idx].append(row.to_dict())
        shard_weights[shard_idx] += float(row["_weight"])

    manifest_rows = []
    for idx, rows in enumerate(shards, start=1):
        shard_name = f"shard_{idx:02d}"
        shard_dir = ensure_dir(outdir / shard_name)
        shard_df = pd.DataFrame(rows).drop(columns=["_weight"], errors="ignore")
        write_tsv(shard_df, shard_dir / "selected_traits.tsv")
        manifest_rows.append(
            {
                "shard_name": shard_name,
                "n_traits": len(shard_df),
                "estimated_weight": int(shard_weights[idx - 1]),
                "selected_path": str(shard_dir / "selected_traits.tsv"),
            }
        )

    write_tsv(pd.DataFrame(manifest_rows), outdir / "shard_manifest.tsv")
    print(f"[done] wrote {len(manifest_rows):,} shard files to {outdir}")


if __name__ == "__main__":
    main()
