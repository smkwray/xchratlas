from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

_shared_data_root = os.environ.get("PROJ_SHARED_DATA_ROOT")
if _shared_data_root:
    DATA_DIR = Path(_shared_data_root).expanduser()
else:
    DATA_DIR = ROOT / "data"

RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
RELEASE_DIR = DATA_DIR / "release"
