"""Shared helpers for the phi-scaling data pipeline.

- HTTP fetch with local cache under data/raw/
- Uniform CSV writer that stamps provenance columns
"""
from __future__ import annotations

import datetime as _dt
import pathlib
import sys
from typing import Optional

import pandas as pd

try:
    import requests
except ImportError:  # pragma: no cover
    requests = None  # type: ignore[assignment]


ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
REFERENCE = ROOT / "data" / "reference"
PROCESSED = ROOT / "data" / "processed"

for _d in (RAW, REFERENCE, PROCESSED):
    _d.mkdir(parents=True, exist_ok=True)


def fetch(url: str, cache_name: str, timeout: int = 30) -> Optional[bytes]:
    """GET `url`, cache raw bytes to data/raw/<cache_name>, return bytes.

    Returns None on any network failure so callers can fall back to the
    shipped reference table.
    """
    cached = RAW / cache_name
    if cached.exists() and cached.stat().st_size > 0:
        return cached.read_bytes()
    if requests is None:
        return None
    try:
        r = requests.get(url, timeout=timeout, headers={"User-Agent": "phi-scaling/1.0"})
        r.raise_for_status()
    except Exception as e:  # noqa: BLE001
        print(f"  [fetch] {url} failed: {e}", file=sys.stderr)
        return None
    cached.write_bytes(r.content)
    return r.content


def write_csv(df: pd.DataFrame, name: str, *, source: str, notes: str = "") -> pathlib.Path:
    """Write df to data/processed/<name> with provenance columns."""
    out = PROCESSED / name
    df = df.copy()
    df["source"] = source
    df["retrieved"] = _dt.date.today().isoformat()
    if notes:
        df["notes"] = notes
    df.to_csv(out, index=False)
    print(f"  wrote {out.relative_to(ROOT)}  ({len(df)} rows)")
    return out
