"""
Test whether spectral line ratios of real elements cluster near phi^k values,
beyond what the hydrogenic 1/n^2 pattern produces.

Approach:
  - For each element, take published wavelengths (NIST ASD, vacuum where noted).
  - Form all pairwise frequency (= energy) ratios r = nu_i / nu_j with r >= 1.
  - Test whether log_phi(r) clusters near integers (phi harmonics).
  - Compare against controls: log_2(r), log_e(r), log_{3/2}(r), log_{sqrt(2)}(r).
  - Null: if no preferred ratio, the fractional part of log_b(r) should be
    uniform on [0, 1). Use Kuiper / Kolmogorov-Smirnov tests.

Interpretation:
  A significant clustering near integers in log_phi but NOT in the other bases
  would be evidence for phi as a preferred harmonic base.
  Clustering in multiple bases (especially simple integer bases) is expected
  from any discrete quantum-mechanical spectrum and is not phi-specific.
"""

import math
import numpy as np
from scipy import stats

PHI = (1 + 5 ** 0.5) / 2  # 1.6180339887...

# -----------------------------------------------------------------------------
# Spectral data (wavelengths in nm, air/vacuum as noted).
# Sources: NIST ASD, CRC Handbook, standard references.
# These are strong, well-documented lines; not an exhaustive list.
# -----------------------------------------------------------------------------

# Hydrogen: Lyman + Balmer + Paschen (air wavelengths; vacuum below 200 nm).
# From NIST ASD.
HYDROGEN = {
    # Lyman (vacuum)
    "Ly-alpha":   121.567,
    "Ly-beta":    102.572,
    "Ly-gamma":    97.253,
    "Ly-delta":    94.974,
    "Ly-epsilon":  93.780,
    # Balmer (air)
    "H-alpha":    656.279,
    "H-beta":     486.135,
    "H-gamma":    434.047,
    "H-delta":    410.174,
    "H-epsilon":  397.007,
    # Paschen (air)
    "Pa-alpha":  1875.10,
    "Pa-beta":   1281.81,
    "Pa-gamma":  1093.81,
}

# Neutral helium (He I) strong lines in air (nm).
HELIUM = {
    "He 388.9":  388.865,
    "He 447.1":  447.148,
    "He 471.3":  471.315,
    "He 492.2":  492.193,
    "He 501.6":  501.568,
    "He 587.6":  587.562,   # D3
    "He 667.8":  667.815,
    "He 706.5":  706.519,
    "He 728.1":  728.135,
    "He 1083":  1083.025,
}

# Neutral sodium (Na I) lines in air (nm). Many from the NIST persistent list.
SODIUM = {
    "Na 330.2":  330.237,
    "Na 330.3":  330.298,
    "Na 568.3":  568.266,
    "Na 568.8":  568.822,
    "Na 589.0 D2": 588.995,
    "Na 589.6 D1": 589.592,
    "Na 615.4":  615.423,
    "Na 616.1":  616.075,
    "Na 819.5":  819.482,
    "Na 819.7":  819.710,
}

# Neutral iron (Fe I) — complex d-shell spectrum.
# Selected strong solar Fraunhofer lines (air, nm).
IRON = {
    "Fe 371.99": 371.9935,
    "Fe 373.49": 373.4864,
    "Fe 374.56": 374.5562,
    "Fe 382.04": 382.0425,
    "Fe 385.99": 385.9911,
    "Fe 404.58": 404.5812,
    "Fe 438.35": 438.3545,
    "Fe 489.07": 489.0754,
    "Fe 495.76": 495.7596,
    "Fe 516.75": 516.7488,
    "Fe 526.95": 526.9537,
    "Fe 532.42": 532.4178,
}

# Neutral neon (Ne I) — strong discharge-tube lines (air, nm).
NEON = {
    "Ne 540.06": 540.0562,
    "Ne 585.25": 585.2488,
    "Ne 588.19": 588.1895,
    "Ne 597.55": 597.5534,
    "Ne 603.00": 603.0000,
    "Ne 607.43": 607.4338,
    "Ne 614.31": 614.3063,
    "Ne 626.65": 626.6495,
    "Ne 633.44": 633.4428,
    "Ne 640.22": 640.2246,
    "Ne 650.65": 650.6528,
    "Ne 659.90": 659.8953,
    "Ne 692.95": 692.9467,
    "Ne 703.24": 703.2413,
}

ELEMENTS = {
    "H":  HYDROGEN,
    "He": HELIUM,
    "Na": SODIUM,
    "Fe": IRON,
    "Ne": NEON,
}

# -----------------------------------------------------------------------------
# Analysis
# -----------------------------------------------------------------------------

def pairwise_ratios(wavelengths_nm):
    """Return array of all pairwise frequency ratios >= 1.
    (Frequencies are 1/lambda; taking the larger over the smaller.)"""
    nu = 1.0 / np.asarray(wavelengths_nm)
    ratios = []
    n = len(nu)
    for i in range(n):
        for j in range(i + 1, n):
            r = nu[i] / nu[j]
            if r < 1.0:
                r = 1.0 / r
            ratios.append(r)
    return np.array(ratios)


def fractional_part(x):
    return x - np.floor(x)


def kuiper_uniform_pvalue(u):
    """Kuiper's test against uniform[0,1). Returns (V, p_approx).
    p approximated with the first few terms of the asymptotic series."""
    u = np.sort(np.asarray(u))
    n = len(u)
    i = np.arange(1, n + 1)
    d_plus = np.max(i / n - u)
    d_minus = np.max(u - (i - 1) / n)
    V = d_plus + d_minus
    # Asymptotic p-value (Stephens 1970)
    lam = (np.sqrt(n) + 0.155 + 0.24 / np.sqrt(n)) * V
    # sum_{k=1..inf} 2*(4*k^2*lam^2 - 1) * exp(-2*k^2*lam^2)
    p = 0.0
    for k in range(1, 101):
        term = 2 * (4 * (k * lam) ** 2 - 1) * math.exp(-2 * (k * lam) ** 2)
        p += term
        if abs(term) < 1e-12:
            break
    p = max(min(p, 1.0), 0.0)
    return V, p


def ks_uniform_pvalue(u):
    """Two-sided KS test of u against Uniform[0,1)."""
    D, p = stats.kstest(u, 'uniform')
    return D, p


def analyze_base(ratios, base, label):
    """Compute fractional parts of log_base(ratios) and test for clustering.
    Returns dict of stats."""
    logs = np.log(ratios) / np.log(base)
    frac = fractional_part(logs)
    # Wrap so that 'near integer' means near 0 OR near 1; measure distance to nearest integer.
    dist_to_int = np.minimum(frac, 1.0 - frac)  # in [0, 0.5]
    ks_D, ks_p = ks_uniform_pvalue(frac)
    kuiper_V, kuiper_p = kuiper_uniform_pvalue(frac)
    # Mean distance to nearest integer; uniform expectation = 0.25.
    mean_d = float(np.mean(dist_to_int))
    # Fraction within 0.05 of an integer; uniform expectation = 0.10.
    frac_near = float(np.mean(dist_to_int < 0.05))
    return {
        "base": base,
        "label": label,
        "n_ratios": len(ratios),
        "mean_dist_to_int": mean_d,
        "uniform_expect_mean": 0.25,
        "frac_within_0.05": frac_near,
        "uniform_expect_frac": 0.10,
        "ks_D": float(ks_D),
        "ks_p": float(ks_p),
        "kuiper_V": float(kuiper_V),
        "kuiper_p": float(kuiper_p),
    }


def analyze_element(name, lines):
    print(f"\n=== {name} ({len(lines)} lines, {len(lines)*(len(lines)-1)//2} pairs) ===")
    labels = list(lines.keys())
    wls = np.array([lines[k] for k in labels])
    ratios = pairwise_ratios(wls)

    bases = [
        (PHI,           "phi   (1.618)"),
        (PHI ** 2,      "phi^2 (2.618)"),
        (2.0,           "2     (octave)"),
        (1.5,           "3/2   (fifth)"),
        (math.e,        "e     (2.718)"),
        (math.sqrt(2),  "sqrt2 (1.414)"),
        (3.0,           "3"),
    ]
    results = []
    for b, lab in bases:
        r = analyze_base(ratios, b, lab)
        results.append(r)

    header = f"  {'base':<16}{'mean d':>9}{'(E=0.25)':>11}{'%<.05':>9}{'(E=.10)':>9}{'KS p':>10}{'Kuiper p':>12}"
    print(header)
    for r in results:
        print(f"  {r['label']:<16}{r['mean_dist_to_int']:>9.4f}{'':>11}"
              f"{r['frac_within_0.05']:>9.3f}{'':>9}"
              f"{r['ks_p']:>10.3g}{r['kuiper_p']:>12.3g}")

    # Also report closest phi^k for top 10 ratios by ... well, we report all.
    logs_phi = np.log(ratios) / np.log(PHI)
    nearest_int = np.round(logs_phi)
    deltas = logs_phi - nearest_int
    # Show distribution in bins of distance-to-integer
    return ratios, results


def summary_across_elements(all_results):
    """Combine KS/Kuiper statistics across elements via Fisher's method."""
    print("\n=== Combined across all elements (Fisher's method on KS p) ===")
    bases_seen = {}
    for elem, results in all_results.items():
        for r in results:
            bases_seen.setdefault(r["label"], []).append(r["ks_p"])
    header = f"  {'base':<16}{'n elements':>12}{'Fisher chi2':>14}{'combined p':>14}"
    print(header)
    for lab, ps in bases_seen.items():
        ps = np.clip(ps, 1e-300, 1.0)
        chi2 = -2.0 * np.sum(np.log(ps))
        dof = 2 * len(ps)
        combined_p = 1.0 - stats.chi2.cdf(chi2, df=dof)
        print(f"  {lab:<16}{len(ps):>12d}{chi2:>14.3f}{combined_p:>14.3g}")


def main():
    print(f"phi = {PHI:.10f}")
    print(f"phi^2 = {PHI**2:.6f}  phi^3 = {PHI**3:.6f}  phi^4 = {PHI**4:.6f}")
    all_results = {}
    for name, lines in ELEMENTS.items():
        _, results = analyze_element(name, lines)
        all_results[name] = results
    summary_across_elements(all_results)

    # Specific phi-harmonic check: count ratios within 1% of phi^k for k=1..6
    print("\n=== Direct count: ratios within 1% of phi^k, k=1..6 ===")
    print("   (Expected by chance ~ 2% per k for uniform log-ratio distribution)")
    print(f"  {'element':<6}{'n pairs':>9}{'k=1':>7}{'k=2':>7}{'k=3':>7}{'k=4':>7}{'k=5':>7}{'k=6':>7}{'any':>7}")
    for name, lines in ELEMENTS.items():
        wls = np.array(list(lines.values()))
        ratios = pairwise_ratios(wls)
        counts = []
        any_mask = np.zeros(len(ratios), dtype=bool)
        for k in range(1, 7):
            target = PHI ** k
            near = np.abs(ratios / target - 1.0) < 0.01
            counts.append(int(near.sum()))
            any_mask |= near
        print(f"  {name:<6}{len(ratios):>9d}" + "".join(f"{c:>7d}" for c in counts)
              + f"{int(any_mask.sum()):>7d}")

    # Same with an integer-base control (base 2 — octaves).
    print("\n=== Control: ratios within 1% of 2^k, k=1..6 ===")
    print(f"  {'element':<6}{'n pairs':>9}{'k=1':>7}{'k=2':>7}{'k=3':>7}{'k=4':>7}{'k=5':>7}{'k=6':>7}{'any':>7}")
    for name, lines in ELEMENTS.items():
        wls = np.array(list(lines.values()))
        ratios = pairwise_ratios(wls)
        counts = []
        any_mask = np.zeros(len(ratios), dtype=bool)
        for k in range(1, 7):
            target = 2.0 ** k
            near = np.abs(ratios / target - 1.0) < 0.01
            counts.append(int(near.sum()))
            any_mask |= near
        print(f"  {name:<6}{len(ratios):>9d}" + "".join(f"{c:>7d}" for c in counts)
              + f"{int(any_mask.sum()):>7d}")


if __name__ == "__main__":
    main()
