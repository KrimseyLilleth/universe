"""
Frequency-ratio analysis across scales.

Question: When we take a collection of natural oscillation frequencies
spanning from atomic to cosmological scales, do the ratios between them
cluster at powers of any particular base -- phi, 2, pi, e -- more than
chance would predict?

Method (circular-statistics approach):

    For a candidate base b, define theta_i = 2*pi * frac(log_b(f_i)).
    Under the null hypothesis that frequencies are drawn log-uniformly
    (no preferred base), the theta_i are uniform on the circle and
    Rayleigh's R = |(1/N) sum_i e^{i theta_i}| has E[R^2] = 1/N.
    The statistic 2 N R^2 is chi-squared(2) under the null, giving
    a clean p-value for "are the log-base-b values clustered mod 1?"

    Kuiper's V is reported alongside Rayleigh because it catches
    multi-modal departures from uniformity that Rayleigh can miss.

    Because "cluster near integer powers of b" is equivalent to
    "fractional parts of log_b(f_i) are concentrated modulo 1,"
    this single test answers the pairwise-ratio question: if every
    f_i is near b^{n_i + c} for integers n_i and a common offset c,
    the circular concentration is maximal.

The frequency list below is a sourced sample of natural oscillations,
not an exhaustive inventory. Any curated list has selection effects,
so we report the test at face value and let the reader judge.
"""

import math
import numpy as np
from scipy import stats

PHI = (1 + 5 ** 0.5) / 2
YEAR = 365.25 * 86400.0  # seconds

# (name, frequency in Hz, category, source note)
FREQUENCIES = [
    # --- cosmological / galactic ---
    ("Hubble frequency (1/age of universe)", 1.0 / (13.8e9 * YEAR), "cosmological",
     "Planck 2018: t_0 ~ 13.8 Gyr"),
    ("Galactic year (Milky Way rotation)", 1.0 / (225e6 * YEAR), "galactic",
     "Sun orbits galactic center in ~225 Myr"),

    # --- solar ---
    ("Solar Schwabe cycle", 1.0 / (11.0 * YEAR), "solar",
     "11-year sunspot cycle"),
    ("Solar p-mode (5 minutes)", 1.0 / 300.0, "solar",
     "Dominant solar oscillation, ~3.3 mHz"),

    # --- planetary / orbital ---
    ("Earth orbital period", 1.0 / YEAR, "planetary", "1 year"),
    ("Lunar sidereal month", 1.0 / (27.3217 * 86400), "planetary",
     "27.32 days"),
    ("Earth rotation (sidereal day)", 1.0 / 86164.1, "planetary",
     "86164.1 s"),

    # --- geophysical (Schumann resonances) ---
    ("Schumann mode 1", 7.83, "geophysical",
     "Fundamental Schumann resonance"),
    ("Schumann mode 2", 14.3, "geophysical", ""),
    ("Schumann mode 3", 20.8, "geophysical", ""),
    ("Schumann mode 4", 27.3, "geophysical", ""),

    # --- biological ---
    ("Circadian rhythm", 1.0 / 86400.0, "biological",
     "Human sleep/wake cycle, ~24h"),
    ("Human breathing (rest)", 0.25, "biological",
     "~15 breaths per minute"),
    ("Human heart rate (rest)", 1.1, "biological",
     "~66 bpm"),
    ("EEG delta (midrange)", 2.5, "biological",
     "Deep sleep, 1-4 Hz"),
    ("EEG theta (midrange)", 6.0, "biological",
     "Drowsy / meditation, 4-8 Hz"),
    ("EEG alpha (midrange)", 10.0, "biological",
     "Relaxed wakefulness, 8-13 Hz"),
    ("EEG beta (midrange)", 21.5, "biological",
     "Active thought, 13-30 Hz"),
    ("EEG gamma (midrange)", 60.0, "biological",
     "Binding / awareness, 30-100 Hz"),

    # --- acoustic ---
    ("Middle C (C4)", 261.6256, "acoustic", "Equal-tempered"),
    ("Concert A (A4)", 440.0, "acoustic", "Standard pitch"),

    # --- atomic / nuclear ---
    ("Hydrogen 21 cm hyperfine", 1.4204058e9, "atomic",
     "Hyperfine transition of neutral hydrogen"),
    ("Cs-133 hyperfine (SI second)", 9.192631770e9, "atomic",
     "SI definition of the second"),
]


def test_base(log_f, base, label):
    """Rayleigh + Kuiper tests on frac(log_base(f)) interpreted as angles."""
    lb = log_f / math.log(base)
    frac = lb - np.floor(lb)            # in [0, 1)
    theta = 2 * math.pi * frac          # in [0, 2 pi)
    N = len(theta)

    # Rayleigh
    C = np.sum(np.cos(theta)) / N
    S = np.sum(np.sin(theta)) / N
    R = math.hypot(C, S)
    rayleigh_stat = 2 * N * R * R
    rayleigh_p = math.exp(-N * R * R) * (
        1 + (2 * N * R * R - (N * R * R) ** 2) / (4 * N)
    )  # small-sample correction; fine for N ~ 20

    # Kuiper's test against uniform on [0, 1)
    kuiper = kuiper_test(frac)

    # Where does the cluster sit?
    mean_angle = math.atan2(S, C)
    if mean_angle < 0:
        mean_angle += 2 * math.pi
    mean_frac = mean_angle / (2 * math.pi)

    return {
        "base": label,
        "R": R,
        "rayleigh_stat": rayleigh_stat,
        "rayleigh_p": rayleigh_p,
        "kuiper_V": kuiper["V"],
        "kuiper_p": kuiper["p"],
        "mean_frac": mean_frac,
    }


def kuiper_test(u):
    """Kuiper's V statistic against U(0,1) with an asymptotic p-value."""
    u = np.sort(np.asarray(u))
    N = len(u)
    i = np.arange(1, N + 1)
    Dplus = np.max(i / N - u)
    Dminus = np.max(u - (i - 1) / N)
    V = (Dplus + Dminus) * (math.sqrt(N) + 0.155 + 0.24 / math.sqrt(N))
    # Stephens (1970) asymptotic p-value
    p = 0.0
    for j in range(1, 101):
        term = (4 * j * j * V * V - 1) * math.exp(-2 * j * j * V * V)
        p += term
    p = 2 * p
    p = max(min(p, 1.0), 0.0)
    return {"V": V, "p": p}


def pairwise_log_ratios(freqs, base):
    """All pairwise log_base(f_i/f_j) for i<j, as distance to nearest integer."""
    lb = np.log(freqs) / math.log(base)
    N = len(lb)
    dists = []
    for i in range(N):
        for j in range(i + 1, N):
            x = lb[i] - lb[j]
            d = abs(x - round(x))
            dists.append(d)
    return np.array(dists)


def main():
    names = [row[0] for row in FREQUENCIES]
    freqs = np.array([row[1] for row in FREQUENCIES], dtype=float)
    log_f = np.log(freqs)
    N = len(freqs)

    print(f"N = {N} frequencies, spanning {freqs.max()/freqs.min():.2e} in range")
    print(f"log10 span: {np.log10(freqs.max()/freqs.min()):.1f} decades\n")

    bases = [
        ("phi  (1.618..)", PHI),
        ("2", 2.0),
        ("e", math.e),
        ("pi", math.pi),
        ("3", 3.0),
        ("10", 10.0),
        ("phi^2 (2.618)", PHI * PHI),
        ("sqrt(2)", math.sqrt(2)),
    ]

    print(f"{'base':20s}  {'R':>6s}  {'Rayleigh p':>11s}  {'Kuiper V':>9s}  {'Kuiper p':>9s}  {'mean frac':>9s}")
    print("-" * 80)
    results = []
    for label, b in bases:
        r = test_base(log_f, b, label)
        results.append(r)
        print(f"{label:20s}  {r['R']:.4f}  {r['rayleigh_p']:11.4f}  "
              f"{r['kuiper_V']:9.4f}  {r['kuiper_p']:9.4f}  {r['mean_frac']:9.4f}")

    # Bonferroni across the 8 bases
    best = min(results, key=lambda x: x["rayleigh_p"])
    print(f"\nBest (smallest Rayleigh p): base {best['base']}, p={best['rayleigh_p']:.4f}")
    print(f"Bonferroni-corrected p (8 bases tested): {min(1.0, best['rayleigh_p']*len(bases)):.4f}")

    # Sanity: the pairwise-ratio view
    print("\nMean distance from nearest integer of log_b(f_i/f_j), all pairs:")
    print(" (uniform null = 0.25; lower = ratios cluster on powers of b)")
    for label, b in bases:
        d = pairwise_log_ratios(freqs, b)
        print(f"  {label:20s}  mean={d.mean():.4f}  median={np.median(d):.4f}")

    # Control: shuffle frequencies log-uniformly over the observed range,
    # repeat many times, see how often phi beats 2.
    print("\nMonte Carlo control (10000 trials of N log-uniform frequencies):")
    rng = np.random.default_rng(0)
    lo, hi = np.log(freqs.min()), np.log(freqs.max())
    phi_wins = 0
    phi_R_null = []
    trials = 10000
    for _ in range(trials):
        sim = np.exp(rng.uniform(lo, hi, size=N))
        r_phi = test_base(np.log(sim), PHI, "phi")["R"]
        r_two = test_base(np.log(sim), 2.0, "2")["R"]
        phi_R_null.append(r_phi)
        if r_phi > r_two:
            phi_wins += 1
    phi_R_null = np.array(phi_R_null)
    observed_R_phi = next(r for r in results if r["base"].startswith("phi "))["R"]
    pct = (phi_R_null >= observed_R_phi).mean()
    print(f"  P(phi R >= observed {observed_R_phi:.4f} | null) = {pct:.4f}")
    print(f"  P(phi beats 2 under null) = {phi_wins/trials:.4f}")

    # Also: drop Schumann higher modes (they are by-construction
    # multiples of the fundamental and inflate any integer-ratio signal).
    keep = [i for i, row in enumerate(FREQUENCIES)
            if row[0] not in ("Schumann mode 2", "Schumann mode 3", "Schumann mode 4")]
    print("\nRepeat test excluding Schumann harmonics 2-4:")
    log_f2 = log_f[keep]
    for label, b in bases:
        r = test_base(log_f2, b, label)
        print(f"  {label:20s}  R={r['R']:.4f}  Rayleigh p={r['rayleigh_p']:.4f}")


if __name__ == "__main__":
    main()
