"""Closed-shell / closed-subshell ionization energies.

Electronic shell closures are local maxima of the first-ionization-energy
curve. The canonical closures are:
  noble gases (n p^6 closure): He, Ne, Ar, Kr, Xe, Rn, Og
  group-12 filled d shell:     Zn, Cd, Hg, Cn
  group-2 filled s shell:      Be, Mg, Ca, Sr, Ba, Ra
  pseudo-closure at N, P, As:  half-filled p^3 (Hund stability bump)

This script derives the table from ionization_energies.csv (produced
by 02_ionization_energies.py) when available, otherwise from the
shipped reference subset. It labels each entry with its closure type.

Sources
-------
Closure assignments: Atkins & de Paula, *Physical Chemistry* 11e, ch.9.
Ionization values: inherited from script 02 (NIST ASD).
"""
from __future__ import annotations

import pandas as pd

from common import PROCESSED, REFERENCE, write_csv

CLOSURES = [
    ("He",  2, "noble_gas"),
    ("Be",  4, "s2"),
    ("N",   7, "half_p3"),
    ("Ne", 10, "noble_gas"),
    ("Mg", 12, "s2"),
    ("P",  15, "half_p3"),
    ("Ar", 18, "noble_gas"),
    ("Ca", 20, "s2"),
    ("As", 33, "half_p3"),
    ("Kr", 36, "noble_gas"),
    ("Sr", 38, "s2"),
    ("Zn", 30, "d10"),
    ("Cd", 48, "d10"),
    ("Sb", 51, "half_p3"),
    ("Xe", 54, "noble_gas"),
    ("Ba", 56, "s2"),
    ("Hg", 80, "d10"),
    ("Bi", 83, "half_p3"),
    ("Rn", 86, "noble_gas"),
    ("Ra", 88, "s2"),
]


def main() -> None:
    print("[03] shell closures")
    ie_path = PROCESSED / "ionization_energies.csv"
    if ie_path.exists():
        ie = pd.read_csv(ie_path)
    else:
        ie = pd.read_csv(REFERENCE / "ionization_energies_reference.csv")

    rows = []
    by_z = {int(r.Z): r for r in ie.itertuples()}
    for sym, Z, kind in CLOSURES:
        if Z not in by_z:
            continue
        r = by_z[Z]
        rows.append({
            "Z": Z,
            "element": getattr(r, "element", sym),
            "closure_type": kind,
            "ionization_energy_eV": float(r.ionization_energy_eV),
        })
    df = pd.DataFrame(rows).sort_values("Z").reset_index(drop=True)
    write_csv(df, "shell_closures.csv",
              source="closure labels: Atkins & de Paula 11e; IEs from NIST ASD",
              notes="IP1 at canonical electronic closures")


if __name__ == "__main__":
    main()
