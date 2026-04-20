"""Nuclear binding energy per nucleon, all known isotopes.

Primary source
--------------
AME2020: Atomic Mass Evaluation 2020
  Wang, Huang, Kondev, Audi, Naimi, Chin. Phys. C 45 (2021) 030003.
  Data file: mass_1.mas20.txt
  Mirror:    https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
  Mirror:    https://amdc.impcas.ac.cn/web/masseval.html

Preprocessing
-------------
The AME2020 mass file is fixed-width. The binding-energy-per-nucleon
column is already tabulated (BIND.ENG/A, in keV). We divide by 1000 to
get MeV and emit one row per (Z, N, A) nuclide.

Reference fallback
------------------
If the live mirror is unreachable, we emit a curated subset of ~60
representative nuclides spanning the chart of nuclides (light, medium,
heavy, superheavy) sourced from the same AME2020 tables as reproduced
in the NNDC NuDat 3.0 database (https://www.nndc.bnl.gov/nudat3/).
"""
from __future__ import annotations

import io
import sys

import pandas as pd

from common import REFERENCE, fetch, write_csv

AME_URL = "https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt"
CACHE = "ame2020_mass.txt"


def _parse_ame(blob: bytes) -> pd.DataFrame:
    text = blob.decode("latin-1", errors="replace")
    rows = []
    for line in text.splitlines():
        # skip header lines (start with '0' flag convention varies);
        # header section ends at a line of dashes. Real data rows start
        # after and have length >= 125. AME2020 columns (see readme):
        #   N, Z, A begin at cols 5, 9, 14;  BIND.ENG/A at cols 55-67 (keV).
        if len(line) < 125 or not line[4:9].strip().isdigit():
            continue
        try:
            N = int(line[4:9])
            Z = int(line[9:14])
            A = int(line[14:19])
            # element symbol, 3 chars
            el = line[20:23].strip()
            be_per_a_str = line[54:67].replace("#", "").replace("*", "").strip()
            if not be_per_a_str:
                continue
            be_per_a_keV = float(be_per_a_str)
        except (ValueError, IndexError):
            continue
        rows.append(
            {
                "Z": Z,
                "N": N,
                "A": A,
                "element": el,
                "binding_energy_per_nucleon_MeV": be_per_a_keV / 1000.0,
            }
        )
    return pd.DataFrame(rows)


def _load_reference() -> pd.DataFrame:
    return pd.read_csv(REFERENCE / "nuclear_binding_reference.csv")


def main() -> None:
    print("[01] nuclear binding energies")
    blob = fetch(AME_URL, CACHE)
    if blob:
        try:
            df = _parse_ame(blob)
            if len(df) < 100:
                raise ValueError(f"only {len(df)} rows parsed, likely format drift")
            source = "AME2020 (Wang et al. 2021), mass_1.mas20.txt"
        except Exception as e:  # noqa: BLE001
            print(f"  parse failed ({e}); using reference table", file=sys.stderr)
            df = _load_reference()
            source = "reference (AME2020 subset via NNDC NuDat 3.0)"
    else:
        print("  no network; using reference table")
        df = _load_reference()
        source = "reference (AME2020 subset via NNDC NuDat 3.0)"

    write_csv(df, "nuclear_binding.csv", source=source,
              notes="BE/A in MeV; higher = more stable")


if __name__ == "__main__":
    main()
