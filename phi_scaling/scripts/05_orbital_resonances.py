"""Orbital-period resonance ratios.

Sources
-------
Solar system periods:
  JPL Planetary Fact Sheet (Williams, D.R., NASA GSFC)
  https://nssdc.gsfc.nasa.gov/planetary/factsheet/
Galilean moons: JPL HORIZONS, Jacobson 2021.
Confirmed exoplanet resonant chains:
  TRAPPIST-1 (Agol et al. 2021, PSJ 2, 1)
  Kepler-223  (Mills et al. 2016, Nature 533, 509)
  Kepler-80   (MacDonald et al. 2016, AJ 152, 105)
  Kepler-90   (Shallue & Vanderburg 2018, AJ 155, 94)
  GJ 876      (Rivera et al. 2010, ApJ 719, 890)
  HR 8799     (Gozdziewski & Migaszewski 2014, MNRAS 440, 3140)
  TOI-178     (Leleu et al. 2021, A&A 649, A26)

Preprocessing
-------------
For each system we list adjacent pairs (P_outer / P_inner) and the
nearest small-integer mean-motion resonance. The processed CSV has one
row per adjacent pair so it can be histogrammed directly.
"""
from __future__ import annotations

import pandas as pd

from common import REFERENCE, write_csv


def main() -> None:
    print("[05] orbital resonance ratios")
    df = pd.read_csv(REFERENCE / "orbital_resonances_reference.csv")
    df["period_ratio"] = df["period_outer_days"] / df["period_inner_days"]
    write_csv(df, "orbital_resonances.csv",
              source="JPL Fact Sheet, JPL HORIZONS, NASA Exoplanet Archive (per-paper citations in CITATIONS.md)",
              notes="period_ratio = P_outer / P_inner for adjacent pairs")


if __name__ == "__main__":
    main()
