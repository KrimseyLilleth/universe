"""Stellar initial-mass-function (IMF) characteristic masses.

The IMF gives dN/dM, the number of stars formed per unit mass. Published
parametrizations have characteristic mass scales: the peak of the
log-normal body, the break between segments of broken power laws, and
the low/high-mass cutoffs. This file compiles those scales across the
major IMF families so they can be tested for phi-scaling.

Sources
-------
Salpeter 1955, ApJ 121, 161                         (single power law)
Miller & Scalo 1979, ApJSS 41, 513                  (log-normal)
Kroupa 2001, MNRAS 322, 231                         (broken power law)
Chabrier 2003, PASP 115, 763                        (log-normal + power law)
Chabrier 2005, ASSL 327, 41                         (system IMF)
Weidner & Kroupa 2006, MNRAS 365, 1333              (IGIMF)
Bastian, Covey & Meyer 2010, ARA&A 48, 339          (review of peak masses)

Entries tagged as 'peak' are log-normal characteristic masses mc
(the median of the underlying log-normal); 'break' entries are power-law
break masses; 'cutoff' entries are the assumed lower/upper ends.
"""
from __future__ import annotations

import pandas as pd

from common import REFERENCE, write_csv


def main() -> None:
    print("[06] stellar IMF characteristic masses")
    df = pd.read_csv(REFERENCE / "stellar_imf_reference.csv")
    write_csv(df, "stellar_imf.csv",
              source="per-row citation (Salpeter 55, MS 79, Kroupa 01, Chabrier 03/05, WK 06)",
              notes="characteristic masses in solar masses")


if __name__ == "__main__":
    main()
