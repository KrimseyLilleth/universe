# Data Source Citations

Every dataset below feeds one CSV in `data/processed/`. Scripts try the
canonical live source first; if the network is unavailable, they fall
back to the curated reference table in `data/reference/`. The `source`
column in each processed CSV records which path was used.

## 1. Nuclear binding energies per nucleon — `nuclear_binding.csv`

- **Primary**: AME2020 (Atomic Mass Evaluation 2020)
  Wang, M., Huang, W.J., Kondev, F.G., Audi, G., Naimi, S.
  *Chinese Physics C* **45**, 030003 (2021).
  File: `mass_1.mas20.txt`
  Mirrors:
  - https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
  - https://amdc.impcas.ac.cn/web/masseval.html
- **Reference fallback**: ~60-nuclide subset reproduced from AME2020
  via NNDC NuDat 3.0 (https://www.nndc.bnl.gov/nudat3/).
- **Preprocessing**: keV → MeV (divide by 1000); one row per (Z, N, A).
  Lines with `#` or `*` flags (systematics / estimated) kept; a `notes`
  column could be added in Phase 2 to flag them.

## 2. First ionization energies — `ionization_energies.csv`

- **Primary**: NIST Atomic Spectra Database (ASD), Ground Levels and
  Ionization Energies.
  Kramida, A., Ralchenko, Yu., Reader, J. and NIST ASD Team (2024).
  https://physics.nist.gov/asd
- **Reference fallback**: Z=1..103 compiled from NIST ASD and
  CRC Handbook of Chemistry and Physics, 104th ed. (2023), §10.
- **Preprocessing**: keep neutral-atom first ionization (spectrum "X I");
  strip uncertainty indicators; eV.

## 3. Electron shell closure energies — `shell_closures.csv`

- Closure assignments: Atkins, P. W. & de Paula, J., *Physical Chemistry*,
  11th ed. (Oxford 2018), ch.9.
- IE values inherited from dataset #2.
- Types: `noble_gas` (np^6), `s2` (ns^2), `d10` (n-1 d^10), `half_p3`.

## 4. Molecular bond dissociation energies — `bond_dissociation.csv`

- **Primary**: Luo, Y.-R., *Comprehensive Handbook of Chemical Bond
  Energies* (CRC Press 2007).
- **Cross-check**: NIST CCCBDB (https://cccbdb.nist.gov/bdelist.asp) and
  NIST Chemistry WebBook (https://webbook.nist.gov/chemistry/).
- **Preprocessing**: homolytic BDE at 298 K in kJ/mol; one row per
  bond type with a representative example molecule.

## 5. Orbital resonance ratios — `orbital_resonances.csv`

- **Solar system periods**: Williams, D. R., NASA/GSFC *Planetary
  Fact Sheet* (https://nssdc.gsfc.nasa.gov/planetary/factsheet/).
- **Galilean moons**: JPL HORIZONS; Jacobson, R. A., satellite ephemerides.
- **Exoplanet chains** (one paper per system):
  - TRAPPIST-1: Agol et al. 2021, *PSJ* **2**, 1.
  - Kepler-223: Mills et al. 2016, *Nature* **533**, 509.
  - Kepler-80: MacDonald et al. 2016, *AJ* **152**, 105.
  - Kepler-90: Shallue & Vanderburg 2018, *AJ* **155**, 94.
  - GJ 876: Rivera et al. 2010, *ApJ* **719**, 890.
  - HR 8799: Goździewski & Migaszewski 2014, *MNRAS* **440**, 3140.
  - TOI-178: Leleu et al. 2021, *A&A* **649**, A26.
- **Preprocessing**: one row per adjacent pair;
  `period_ratio = P_outer / P_inner` computed at load time.

## 6. Stellar IMF characteristic masses — `stellar_imf.csv`

- Salpeter, E. E. 1955, *ApJ* **121**, 161.
- Miller, G. E. & Scalo, J. M. 1979, *ApJSS* **41**, 513.
- Kroupa, P. 2001, *MNRAS* **322**, 231.
- Chabrier, G. 2003, *PASP* **115**, 763.
- Chabrier, G. 2005, *ASSL* **327**, 41.
- Weidner, C. & Kroupa, P. 2006, *MNRAS* **365**, 1333.
- Bastian, N., Covey, K. R. & Meyer, M. R. 2010, *ARA&A* **48**, 339.
- **Preprocessing**: one row per published scale with `scale_type` in
  {peak, break, cutoff_low, cutoff_high}; mass in solar masses.

## 7. Biological scales — `biological_scales.csv`

Long format; `modality` column is one of `cell_size`, `neural_band`,
`hrv_band`.

- **Cell sizes**:
  - BioNumbers: Milo, R. et al. 2010, *NAR* **38**, D750.
    https://bionumbers.hms.harvard.edu/
  - Milo, R. & Phillips, R., *Cell Biology by the Numbers*
    (Garland 2015).
- **Neural oscillation bands**:
  - Buzsáki, G. & Draguhn, A. 2004, *Science* **304**, 1926.
  - Niedermeyer, E. & Lopes da Silva, F., *Electroencephalography:
    Basic Principles, Clinical Applications, and Related Fields*,
    6th ed. (Lippincott 2011).
- **Heart-rate-variability bands**:
  - Task Force of the ESC/NASPE 1996, *Circulation* **93**, 1043.
  - Shaffer, F. & Ginsberg, J. P. 2017, *Front. Public Health* **5**, 258.
- **Preprocessing**: `value` is the band-center / typical value; where
  applicable `value_low` and `value_high` give the band edges.
