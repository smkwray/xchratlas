from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

from .io import choose_column
from .region import classify_x_region

GWS_P = 5e-8
GWS_NEGLOG10 = -math.log10(GWS_P)


PVAL_PRIORITY = [
    "neglog10_pval_meta_hq",
    "neglog10_pval_meta",
    "pval_meta_hq",
    "pval_meta",
    "log_pval_meta_hq",
    "log_pval_meta",
]


def choose_best_metric_column(columns: list[str], preferred: list[str]) -> str | None:
    return choose_column(columns, exact=preferred)


def choose_best_pvalue_column(columns: list[str]) -> str | None:
    hit = choose_column(columns, exact=PVAL_PRIORITY)
    if hit:
        return hit
    candidates = [c for c in columns if "pval" in c.lower()]
    if not candidates:
        candidates = [c for c in columns if "pvalue" in c.lower()]
    if not candidates:
        candidates = [c for c in columns if "log" in c.lower() and "p" in c.lower()]
    if not candidates:
        return None
    # Prefer meta + hq, then meta, then anything.
    candidates = sorted(
        candidates,
        key=lambda c: (
            "meta" not in c.lower(),
            "hq" not in c.lower(),
            c.lower(),
        ),
    )
    return candidates[0]


def infer_pvalue_scale(series: pd.Series, column_name: str) -> str:
    col = column_name.lower()
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        raise ValueError(f"Could not infer p-value scale for empty series: {column_name}")

    if "neglog10" in col:
        return "neglog10"
    if "ln" in col:
        return "ln"
    if s.max() <= 1.0 and s.min() >= 0.0:
        return "raw"
    if "log" in col and s.min() >= 0:
        return "neglog10"
    if s.max() <= 0:
        return "ln"
    # Fall back to raw if values look bounded, else neglog10 if large positive.
    if s.max() <= 1.0:
        return "raw"
    return "neglog10"


def to_raw_pvalue(series: pd.Series, column_name: str) -> tuple[pd.Series, str]:
    numeric = pd.to_numeric(series, errors="coerce")
    scale = infer_pvalue_scale(numeric, column_name)
    if scale == "raw":
        p = numeric
    elif scale == "neglog10":
        p = 10 ** (-numeric)
    elif scale == "ln":
        p = np.exp(numeric)
    else:
        raise ValueError(f"Unsupported p-value scale: {scale}")

    p = p.clip(lower=1e-300, upper=1.0)
    return p, scale


def to_neglog10_pvalue(series: pd.Series, column_name: str) -> tuple[pd.Series, str]:
    p, scale = to_raw_pvalue(series, column_name)
    return -np.log10(p), scale


def _pick_optional(df: pd.DataFrame, choices: list[str]) -> str | None:
    for choice in choices:
        if choice in df.columns:
            return choice
    return None


def call_loci_from_dataframe(
    df: pd.DataFrame,
    *,
    trait_id: str,
    window_bp: int = 500_000,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if "pos" not in df.columns:
        raise ValueError("DataFrame must include a 'pos' column.")
    if "chr" not in df.columns:
        raise ValueError("DataFrame must include a 'chr' column.")

    pcol = choose_best_pvalue_column(list(df.columns))
    if pcol is None:
        raise ValueError(f"Could not find a usable p-value column for trait {trait_id}")

    d = df.copy()
    d["pos"] = pd.to_numeric(d["pos"], errors="coerce")
    d = d.dropna(subset=["pos"]).copy()
    d["pos"] = d["pos"].astype(int)
    d["x_region"] = d["pos"].map(classify_x_region)

    raw_p, scale = to_raw_pvalue(d[pcol], pcol)
    d["lead_pvalue"] = raw_p
    d["lead_neglog10_pvalue"] = -np.log10(raw_p)

    beta_col = _pick_optional(d, ["beta_meta_hq", "beta_meta", "beta"])
    se_col = _pick_optional(d, ["se_meta_hq", "se_meta", "se"])

    sig = d.loc[d["lead_pvalue"] < GWS_P].copy()
    sig = sig.sort_values(["x_region", "pos", "lead_neglog10_pvalue"], ascending=[True, True, False])

    if sig.empty:
        summary = {
            "trait_id": trait_id,
            "n_sig_variants": 0,
            "n_loci": 0,
            "lead_neglog10_pvalue": np.nan,
            "pvalue_column": pcol,
            "pvalue_scale": scale,
            "par_only_signal": False,
            "any_nonpar_locus": False,
        }
        return pd.DataFrame(), summary

    loci: list[dict[str, Any]] = []
    current_rows: list[pd.Series] = []
    current_region: str | None = None
    current_start: int | None = None
    current_end: int | None = None

    def flush_current() -> None:
        nonlocal current_rows, current_region, current_start, current_end
        if not current_rows:
            return
        block = pd.DataFrame(current_rows)
        lead = block.sort_values("lead_neglog10_pvalue", ascending=False).iloc[0]
        locus_id = f"{trait_id}:{current_region}:{int(current_start)}-{int(current_end)}"
        row = {
            "trait_id": trait_id,
            "locus_id": locus_id,
            "x_region": current_region,
            "locus_start": int(current_start),
            "locus_end": int(current_end),
            "n_sig_variants": int(len(block)),
            "lead_chr": str(lead.get("chr", "X")),
            "lead_pos": int(lead["pos"]),
            "lead_ref": lead.get("ref"),
            "lead_alt": lead.get("alt"),
            "lead_rsid": lead.get("rsid"),
            "lead_varid": lead.get("varid"),
            "lead_beta": lead.get(beta_col) if beta_col else np.nan,
            "lead_se": lead.get(se_col) if se_col else np.nan,
            "lead_pvalue": float(lead["lead_pvalue"]),
            "lead_neglog10_pvalue": float(lead["lead_neglog10_pvalue"]),
            "pvalue_column": pcol,
            "pvalue_scale": scale,
        }
        loci.append(row)
        current_rows = []
        current_region = None
        current_start = None
        current_end = None

    for _, row in sig.iterrows():
        pos = int(row["pos"])
        region = row["x_region"]
        if not current_rows:
            current_rows = [row]
            current_region = region
            current_start = pos
            current_end = pos
            continue

        assert current_end is not None
        if region == current_region and pos - current_end <= window_bp:
            current_rows.append(row)
            current_end = pos
        else:
            flush_current()
            current_rows = [row]
            current_region = region
            current_start = pos
            current_end = pos

    flush_current()
    loci_df = pd.DataFrame(loci).sort_values(["x_region", "lead_pos"]).reset_index(drop=True)

    summary = {
        "trait_id": trait_id,
        "n_sig_variants": int(len(sig)),
        "n_loci": int(len(loci_df)),
        "lead_neglog10_pvalue": float(loci_df["lead_neglog10_pvalue"].max()),
        "pvalue_column": pcol,
        "pvalue_scale": scale,
        "par_only_signal": bool(set(loci_df["x_region"]) <= {"PAR1", "PAR2"}),
        "any_nonpar_locus": bool((loci_df["x_region"] == "nonPAR").any()),
    }
    return loci_df, summary
