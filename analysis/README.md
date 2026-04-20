# Spectral φ-Harmonic Test

## Hypothesis

If the golden ratio φ ≈ 1.618 governs atomic shell geometry, then the
ratios between spectral lines of a given element should show φ-related
harmonic relationships — ratios clustering near φᵏ — distinct from
(or in addition to) the hydrogen-series 1/n² pattern.

## Data

Published NIST Atomic Spectra Database wavelengths for:

| Element | Lines | Pairwise ratios |
|---|---|---|
| H (Lyman + Balmer + Paschen) | 13 | 78 |
| He I  | 10 | 45 |
| Na I  | 10 | 45 |
| Fe I  | 12 | 66 |
| Ne I  | 14 | 91 |

All wavelengths hard-coded in `spectral_phi_analysis.py` with sources noted.

## Three tests

### Test 1 — Fixed-base clustering (`spectral_phi_analysis.py`)

For each pairwise frequency ratio r ≥ 1, compute the fractional part of
log_b(r) for several candidate bases. If b is a privileged harmonic
base, the fractional parts should cluster near 0 (ratios near bᵏ).
Kuiper and Kolmogorov–Smirnov tests against Uniform[0,1).

**Result:** Every element rejects uniformity for nearly every base tested
(φ, φ², 2, e, √2, 3, 3/2). This is *expected*: any discrete quantum
spectrum gives structured, non-random ratios. φ is not singled out.

### Test 2 — Base scan (`spectral_phi_scan.py`)

Scan candidate bases b from 1.05 to 3.00 in 0.005 steps. For each b,
compute the Rayleigh Z-statistic for concentration of fractional
log_b(r) near integers. If φ is privileged, it should be a local peak
of Z. Result across all elements combined:

| Base | Combined Rayleigh sum-Z | Rank (of 391 grid points) |
|---|---|---|
| φ  (1.618)  | 68.9  | **258 / 391** |
| φ² (2.618)  | 121.2 |  77 / 391 |
| 2           |  92.1 | 201 / 391 |
| e           | 126.2 |  57 / 391 |
| 3/2         |  48.9 | 324 / 391 |
| √2          |  48.8 | 325 / 391 |

φ ranks in the **bottom third** — it is not preferred. The apparent
strength of large bases (e, 3) is an artifact of the finite wavelength
window: when the ratios lie in a narrow range, log_b(r) for large b
automatically piles up near zero. Hydrogen alone shows a weak peak at
b ≈ 1.615 (Z = 9.09 vs. 8.99 at φ), but this is within Monte-Carlo
noise once the ~400-point multiple-comparison correction is applied.

### Test 3 — Direct shell-radius prediction (`phi_shell_vs_bohr.py`)

The most concrete reading of "φ governs shell geometry" is that
shell radii scale as rₙ = r₀·φⁿ, giving energy levels Eₙ ∝ −1/φⁿ
instead of the Bohr/Schrödinger Eₙ ∝ −1/n². Both models have one free
constant; anchor each to Lyman-α and predict Ly-β through Ly-ε:

| n | Observed λ (nm) | 1/n² pred | err % | 1/φⁿ pred | err % |
|---|---|---|---|---|---|
| 2 | 121.567 | 121.567 | anchor | 121.567 | anchor |
| 3 | 102.572 | 102.572 | +0.000 |  75.13 | −26.75 |
| 4 |  97.253 |  97.254 | +0.001 |  60.78 | −37.50 |
| 5 |  94.974 |  94.974 | +0.000 |  54.37 | −42.76 |
| 6 |  93.780 |  93.780 | +0.000 |  51.04 | −45.58 |

RMS relative error (excluding the anchor):

- **Bohr 1/n² model: 0.0004 %**
- **φ-shell 1/φⁿ model: 38.8 %**
- φ-shell is ~10⁵× worse.

## Conclusion

The spectral-line data reject the literal reading of "φ-spacing
governs shell geometry." Specifically:

1. Pairwise line ratios across H, He, Na, Fe, Ne do not cluster near
   φᵏ more strongly than near φ-neighboring bases (φ is 258ᵗʰ of 391).
2. Hydrogen energy levels follow 1/n² to 0.0004 % RMS; a φⁿ shell
   spacing misses by ~40 %.
3. Any non-uniformity found in log-ratio fractional parts is equally
   present for integer and transcendental bases — a generic consequence
   of discrete quantum spectra, not evidence for φ.

The torus-geometry paper's core result — D = T·φ⁴ at cosmological
scales — stands or falls on its own data. It does **not** appear to
generalize into a φⁿ shell structure inside atoms, and claims in that
direction are not supported by this dataset.

## Reproducing

```bash
cd analysis
python3 spectral_phi_analysis.py
python3 spectral_phi_scan.py
python3 phi_shell_vs_bohr.py
```

Dependencies: `numpy`, `scipy`.
