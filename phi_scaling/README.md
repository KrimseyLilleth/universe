# Phi-Scaling in Natural Stability Thresholds

Testing a hypothesis from the Unified Physics framework: that stable
configurations across physical systems cluster at scales related by the
golden ratio φ ≈ 1.618034 or its integer powers.

This subproject is **Phase 1: data compilation**. It pulls seven datasets
from canonical public sources and writes cleaned CSVs to `data/processed/`.
Every CSV includes provenance columns (`source`, `retrieved`) and every
script documents its citation at the top.

## Layout

```
phi_scaling/
├── README.md                     this file
├── requirements.txt              Python deps (pandas, numpy, requests)
├── data/
│   ├── raw/                      raw downloads cached here
│   ├── reference/                literature-compiled fallbacks (versioned)
│   ├── processed/                clean CSVs produced by scripts (gitignored)
│   └── CITATIONS.md              full source list + preprocessing notes
└── scripts/
    ├── common.py                 shared helpers (cache, CSV writer)
    ├── 01_nuclear_binding.py     AME2020 binding energy per nucleon
    ├── 02_ionization_energies.py NIST ASD first ionization energies
    ├── 03_shell_closures.py      noble-gas / closed-shell ionization
    ├── 04_bond_dissociation.py   NIST webbook BDEs
    ├── 05_orbital_resonances.py  solar system + exoplanet period ratios
    ├── 06_stellar_imf.py         IMF characteristic masses
    ├── 07_biological_scales.py   cell sizes, EEG bands, HRV bands
    └── run_all.py                orchestrator
```

## Running

```bash
python3 -m pip install -r requirements.txt
python3 scripts/run_all.py
```

Each script is independent and can be run standalone. Scripts first attempt
a live fetch from the canonical source into `data/raw/`; if the network is
unavailable they fall back to the curated tables in `data/reference/`.
The output CSV records which path was used in its `source` column so
downstream analysis can filter on provenance.

## Phase 2 (not in this commit)

Binned-log-histogram and KDE peak detection, with a null model that
randomizes log-scale positions, to test whether peak separations cluster
near ln(φ) ≈ 0.481.
