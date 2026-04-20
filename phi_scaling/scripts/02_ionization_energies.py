"""First ionization energy (IP1) for every element Z=1..118.

Primary source
--------------
NIST Atomic Spectra Database (ASD), Ground Levels and Ionization Energies.
  Kramida, A., Ralchenko, Yu., Reader, J. and NIST ASD Team (2024).
  https://physics.nist.gov/asd  (specifically the 'Ionization Energies
  Data' page, form submission with spectrum=All).
  Endpoint used below:
  https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=H-Og&units=1&format=2&output=0

Preprocessing
-------------
The ASD CGI returns CSV when format=2; IE is in eV. We keep only the
neutral atom IE (i.e. the first ionization energy), discarding
spectrum strings past the first ionization stage. Uncertainties and
asterisks are stripped.

Reference fallback
------------------
Hand-compiled from the same NIST ASD (CODATA-recommended rounded
values; see also CRC Handbook 104th ed., §10).
"""
from __future__ import annotations

import io
import re
import sys

import pandas as pd

from common import REFERENCE, fetch, write_csv

NIST_URL = (
    "https://physics.nist.gov/cgi-bin/ASD/ie.pl"
    "?spectra=H-Og&units=1&format=2&output=0&at_num_out=on&el_name_out=on"
    "&seq_out=on&shells_out=on&level_out=on&e_out=0&submit=Retrieve+Data"
)
CACHE = "nist_asd_ie.csv"


def _parse_nist(blob: bytes) -> pd.DataFrame:
    text = blob.decode("utf-8", errors="replace")
    if "Input Error" in text or "<html" in text.lower()[:500] and "<pre" not in text.lower():
        raise ValueError("NIST returned an error page, not CSV")
    # NIST wraps CSV in <pre>...</pre>; strip tags
    text = re.sub(r"<[^>]+>", "", text)
    df = pd.read_csv(io.StringIO(text))
    df.columns = [c.strip().strip('"') for c in df.columns]
    want = {"At. num": "Z", "El. name": "element", "Sp. Name": "spectrum",
            "Ionization Energy (eV)": "ionization_energy_eV"}
    df = df[[c for c in want if c in df.columns]].rename(columns=want)
    # keep only neutral (first) ionization: spectrum "X I"
    df = df[df["spectrum"].astype(str).str.match(r"^\s*\w+\s+I\s*$")]
    df["ionization_energy_eV"] = (
        df["ionization_energy_eV"].astype(str)
        .str.replace(r"[()\[\]a-zA-Z*]", "", regex=True)
        .str.strip()
        .astype(float)
    )
    return df[["Z", "element", "ionization_energy_eV"]].sort_values("Z").reset_index(drop=True)


def _load_reference() -> pd.DataFrame:
    return pd.read_csv(REFERENCE / "ionization_energies_reference.csv")


def main() -> None:
    print("[02] first ionization energies")
    blob = fetch(NIST_URL, CACHE)
    if blob:
        try:
            df = _parse_nist(blob)
            if len(df) < 80:
                raise ValueError(f"only {len(df)} rows parsed")
            source = "NIST ASD (Kramida et al. 2024)"
        except Exception as e:  # noqa: BLE001
            print(f"  parse failed ({e}); using reference", file=sys.stderr)
            df = _load_reference()
            source = "reference (NIST ASD, CRC Handbook 104th ed.)"
    else:
        print("  no network; using reference table")
        df = _load_reference()
        source = "reference (NIST ASD, CRC Handbook 104th ed.)"
    write_csv(df, "ionization_energies.csv", source=source,
              notes="IP1 of neutral atom, ground state, in eV")


if __name__ == "__main__":
    main()
