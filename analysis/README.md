# Spectral φ-Harmonic Test

## Question

Do pairwise ratios between spectral lines of real elements cluster
near φᵏ, as would be expected if φ ≈ 1.618 were a preferred harmonic
base in atomic structure?

## Data

Published NIST Atomic Spectra Database wavelengths for five elements:

| Element | Lines | Pairwise ratios |
|---|---|---|
| H (Lyman + Balmer + Paschen) | 13 | 78 |
| He I  | 10 | 45 |
| Na I  | 10 | 45 |
| Fe I  | 12 | 66 |
| Ne I  | 14 | 91 |

Hard-coded in `spectral_phi_analysis.py` with sources noted. No model
is imposed — only the published wavelengths are used.

## Tests (all three are just ways of looking at the same ratios)

### 1. Direct φᵏ hit counts (`phi_ratio_hits.py`)

For each of 325 pairwise frequency ratios, check whether it falls
within tolerance of φᵏ for k ∈ {−4…+4}\{0}. Compare to an empirical
null that preserves the ratio *distribution* but randomly shifts
each log-ratio by a uniform offset in [−log φ/2, +log φ/2] (500
trials). This removes any quantum structure without changing the
overall density of ratios.

| Tolerance | Observed hits | Empirical null (mean ± σ) | z-score |
|---|---|---|---|
| 0.5 % |  4 |  2.6 ± 1.6 | +0.88 |
| 1.0 % |  7 |  5.4 ± 2.3 | +0.69 |
| 2.0 % | 13 | 10.8 ± 3.1 | +0.72 |

**No significant excess.** Observed counts are within ~1σ of the
null expectation.

### 2. Per-element breakdown

Every hit at 1 % tolerance:

```
H:  Ly-ε / H-ε       ratio = 4.2334   phi^+3 (err +0.063%)
H:  Ly-β / H-γ       ratio = 4.2316   phi^+3 (err +0.105%)
H:  Ly-γ / H-δ       ratio = 4.2176   phi^+3 (err +0.436%)
H:  H-β  / Pa-β      ratio = 2.6367   phi^+2 (err +0.714%)
H:  Ly-δ / H-α       ratio = 6.9101   phi^+4 (err +0.817%)
He: 667.8 / 1083     ratio = 1.6217   phi^+1 (err +0.229%)
He: 447.1 / 728.1    ratio = 1.6284   phi^+1 (err +0.641%)
Na: (none)
Fe: (none)
Ne: (none)
```

Three of the five elements produce zero hits. The Hydrogen hits cluster
near φ³ ≈ 4.236, but this is a coincidence of the Rydberg formula:
for large n,m the Lyman-to-Balmer frequency ratio

(1 − 1/n²) / (¼ − 1/m²)  →  4  as n,m → ∞.

φ³ = 4.236 happens to sit 5.8 % above 4, so many (n,m) pairs land near
it without φ playing any role. For example 800/189, 6860/1620, 2160/512
are the exact rationals underlying the top three H hits.

### 3. Continuous base scan (`spectral_phi_scan.py`)

If φ were privileged, the clustering statistic (Rayleigh Z of fractional
log_b(ratio)) should peak at b = φ when we scan b from 1.05 to 3.00.
Combined across all five elements:

| Base | Rayleigh sum-Z | Rank (of 391) |
|---|---|---|
| φ  (1.618)  |  68.9 | **258** |
| φ² (2.618)  | 121.2 |  77 |
| 2           |  92.1 | 201 |
| e           | 126.2 |  57 |
| 3/2         |  48.9 | 324 |
| √2          |  48.8 | 325 |

φ sits in the bottom third of candidate bases. The apparent preference
for large bases (e, 3) is a finite-window artifact: when ratios live in
a narrow range, log_b with large b compresses them near zero.

## Conclusion

Ratios in existing spectral data do **not** show φ-related harmonic
relationships beyond what you'd get by chance. The ~7 ostensibly
φ-close hits at 1 % tolerance are consistent with an empirical null
(z = +0.69), and the ones that do occur are traceable to rational
accidents of the Rydberg formula rather than any φ-structure. Three
of the five elements give zero hits.

## Reproducing

```bash
cd analysis
python3 spectral_phi_analysis.py   # fixed-base clustering tests
python3 spectral_phi_scan.py       # continuous base scan
python3 phi_ratio_hits.py          # direct phi^k hit counts vs null
```

Dependencies: `numpy`, `scipy`.
