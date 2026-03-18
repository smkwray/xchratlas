from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd


def ensure_dir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def read_tsv(path_or_buf, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path_or_buf, sep="\t", **kwargs)


def write_tsv(df: pd.DataFrame, path: str | Path, index: bool = False) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    df.to_csv(path, sep="\t", index=index, compression=compression)


def coerce_bool(value) -> bool | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    s = str(value).strip().lower()
    if s in {"true", "1", "yes", "y", "t"}:
        return True
    if s in {"false", "0", "no", "n", "f"}:
        return False
    return None


def choose_column(
    columns: Sequence[str],
    *,
    exact: Sequence[str] | None = None,
    contains_all: Sequence[str] | None = None,
    contains_any: Sequence[str] | None = None,
) -> str | None:
    lowered = {c.lower(): c for c in columns}
    if exact:
        for item in exact:
            if item.lower() in lowered:
                return lowered[item.lower()]

    if contains_all:
        for c in columns:
            cl = c.lower()
            if all(token.lower() in cl for token in contains_all):
                return c

    if contains_any:
        for c in columns:
            cl = c.lower()
            if any(token.lower() in cl for token in contains_any):
                return c
    return None


def open_text_auto(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def build_search_blob(df: pd.DataFrame, candidate_columns: Sequence[str]) -> pd.Series:
    usable = [c for c in candidate_columns if c in df.columns]
    if not usable:
        usable = list(df.select_dtypes(include=["object"]).columns)
    if not usable:
        return pd.Series([""] * len(df), index=df.index)
    blob = df[usable].fillna("").astype(str).agg(" | ".join, axis=1)
    return blob.str.lower()


def maybe_read_bgz_tsv(path: str | Path, nrows: int | None = None) -> pd.DataFrame:
    path = Path(path)
    compression = "gzip" if path.suffix.lower() in {".gz", ".bgz"} else "infer"
    return pd.read_csv(path, sep="\t", compression=compression, nrows=nrows)


def text_lines_from_bytes_stream(byte_stream: Iterable[bytes]) -> Iterable[str]:
    for chunk in byte_stream:
        if not chunk:
            continue
        yield chunk.decode("utf-8")
