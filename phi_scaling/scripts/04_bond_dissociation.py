"""Bond dissociation energies (BDE) for common covalent bonds.

Primary source
--------------
NIST Chemistry WebBook, Computational Chemistry Comparison and Benchmark
DataBase (CCCBDB) and Luo, *Comprehensive Handbook of Chemical Bond
Energies* (CRC Press 2007), values cross-referenced with:
  https://webbook.nist.gov/chemistry/
  https://cccbdb.nist.gov/bdelist.asp

Preprocessing
-------------
BDE in kJ/mol at 298 K (standard enthalpy of homolytic dissociation).
The WebBook does not expose a single bulk CSV for diatomic BDEs, so we
ship a curated reference table covering the most common bond types
(single, double, triple; C-H, C-C, C-O, C-N, O-H, N-H, etc.) and
selected diatomics (H2, N2, O2, F2, CO, NO). The table is derived from
Luo 2007 (primary) with WebBook cross-check.

This script currently emits the reference table directly; a live-fetch
path would need per-bond queries against cccbdb.nist.gov.
"""
from __future__ import annotations

import pandas as pd

from common import REFERENCE, write_csv


def main() -> None:
    print("[04] bond dissociation energies")
    df = pd.read_csv(REFERENCE / "bond_dissociation_reference.csv")
    write_csv(df, "bond_dissociation.csv",
              source="Luo 2007 (CRC); cross-checked against NIST CCCBDB",
              notes="homolytic BDE at 298 K, kJ/mol")


if __name__ == "__main__":
    main()
