"""Biological scale data across three modalities.

Three sub-datasets are produced, then concatenated in a long-format CSV
with a `modality` column so downstream analysis can slice on it:

  modality = cell_size       (characteristic diameters, micrometers)
  modality = neural_band     (EEG/LFP rhythm frequency bands, Hz)
  modality = hrv_band        (heart-rate-variability spectral bands, Hz)

Sources
-------
Cell sizes:
  BioNumbers database (Milo et al. 2010, NAR 38 D750)
    https://bionumbers.hms.harvard.edu/
  Milo & Phillips, *Cell Biology by the Numbers* (Garland 2015).

Neural oscillation bands:
  Buzsaki & Draguhn 2004, Science 304, 1926 (canonical band table)
  Niedermeyer & Lopes da Silva, *Electroencephalography* 6e (2011).

Heart-rate-variability bands:
  Task Force of the European Society of Cardiology and the North American
  Society of Pacing and Electrophysiology, Circulation 93 (1996) 1043.
  Shaffer & Ginsberg 2017, Front. Public Health 5, 258.
"""
from __future__ import annotations

import pandas as pd

from common import REFERENCE, write_csv


def main() -> None:
    print("[07] biological scales")
    df = pd.read_csv(REFERENCE / "biological_scales_reference.csv")
    write_csv(df, "biological_scales.csv",
              source="BioNumbers / Milo&Phillips 2015; Buzsaki&Draguhn 2004; Task Force 1996 / Shaffer&Ginsberg 2017",
              notes="long format with `modality` column; see unit column")


if __name__ == "__main__":
    main()
