"""
Sharper test: scan a continuous grid of candidate bases b in [1.05, 3.0]
and ask whether phi = 1.618 sits at a local extremum of the clustering
statistic. If phi is privileged, it should be a LOCAL MAXIMUM of
clustering-near-integers.

For each base b we compute:
  fractional parts of log_b(r) for all pairwise frequency ratios r (>=1),
  then the Rayleigh statistic Z = |<exp(i*2*pi*frac)>|^2 * n.

Large Z => strong concentration on the unit circle (clustering near an
integer phase). A privileged base appears as a peak in Z vs b.
"""

import math
import numpy as np
from scipy import stats

PHI = (1 + 5 ** 0.5) / 2

# Same dataset as the main script (imported inline to keep this standalone).
from spectral_phi_analysis import ELEMENTS, pairwise_ratios


def rayleigh_Z(ratios, base):
    logs = np.log(ratios) / np.log(base)
    frac = logs - np.floor(logs)
    theta = 2 * np.pi * frac
    c = np.mean(np.cos(theta))
    s = np.mean(np.sin(theta))
    R2 = c * c + s * s
    n = len(frac)
    # Rayleigh Z = n*R^2 ; p ~ exp(-Z) for large n.
    Z = n * R2
    p = math.exp(-Z) if Z < 700 else 0.0
    return Z, p


def scan_base(ratios, bases):
    Zs = np.zeros_like(bases)
    ps = np.zeros_like(bases)
    for i, b in enumerate(bases):
        Z, p = rayleigh_Z(ratios, b)
        Zs[i] = Z
        ps[i] = p
    return Zs, ps


def local_max_report(bases, Zs, target, window=0.15):
    """Check whether `target` is a local maximum of Zs within +/- window."""
    mask = np.abs(bases - target) <= window
    if not mask.any():
        return None
    local_Z_max = Zs[mask].max()
    # Find base that achieves the max within the window.
    best_b = bases[mask][np.argmax(Zs[mask])]
    # Z at target (closest grid point).
    idx_target = int(np.argmin(np.abs(bases - target)))
    Z_at_target = Zs[idx_target]
    return {
        "target": float(target),
        "Z_at_target": float(Z_at_target),
        "best_base_in_window": float(best_b),
        "best_Z_in_window": float(local_Z_max),
        "is_local_peak": bool(abs(best_b - target) < (bases[1] - bases[0]) * 2),
    }


def main():
    bases = np.linspace(1.05, 3.00, 391)  # step ~0.005

    print("Scanning base b from 1.05 to 3.00 in steps of ~0.005.")
    print("Rayleigh Z measures clustering of fractional log_b(ratio) near integers.")
    print(f"phi = {PHI:.6f}\n")

    combined_Zs = np.zeros_like(bases)
    per_elem = {}
    for name, lines in ELEMENTS.items():
        wls = np.array(list(lines.values()))
        ratios = pairwise_ratios(wls)
        Zs, ps = scan_base(ratios, bases)
        # Fisher-combine across elements: equivalent to summing Z's (since p ~ exp(-Z))
        combined_Zs += Zs
        per_elem[name] = Zs

        rep = local_max_report(bases, Zs, PHI, window=0.15)
        global_best_idx = int(np.argmax(Zs))
        print(f"{name}: global peak at b = {bases[global_best_idx]:.3f} (Z={Zs[global_best_idx]:.2f});"
              f" at phi Z={rep['Z_at_target']:.2f};"
              f" best near phi (±0.15) at b={rep['best_base_in_window']:.3f}"
              f" (Z={rep['best_Z_in_window']:.2f})")

    print()
    rep = local_max_report(bases, combined_Zs, PHI, window=0.15)
    idx = int(np.argmax(combined_Zs))
    print(f"COMBINED: global peak at b = {bases[idx]:.3f} (sumZ={combined_Zs[idx]:.2f})")
    print(f"COMBINED: at phi sumZ = {rep['Z_at_target']:.2f}")
    print(f"COMBINED: best base within ±0.15 of phi: b = {rep['best_base_in_window']:.3f}"
          f" (sumZ = {rep['best_Z_in_window']:.2f})")

    # Top 10 peaks in combined signal.
    order = np.argsort(-combined_Zs)
    print("\nTop 15 combined peaks:")
    shown = 0
    seen = []
    for i in order:
        b = bases[i]
        if any(abs(b - s) < 0.03 for s in seen):
            continue
        seen.append(b)
        print(f"  b = {b:.3f}   sumZ = {combined_Zs[i]:.2f}")
        shown += 1
        if shown >= 15:
            break

    # Ranks of phi, phi^2, 2, e, 3 in the combined scan.
    print("\nRank (lower = stronger) of specific bases in the combined scan:")
    for name, target in [("phi", PHI), ("phi^2", PHI**2), ("e", math.e),
                         ("2", 2.0), ("3/2", 1.5), ("sqrt(2)", math.sqrt(2))]:
        idx_t = int(np.argmin(np.abs(bases - target)))
        Z_t = combined_Zs[idx_t]
        rank = int(np.sum(combined_Zs > Z_t)) + 1
        print(f"  {name:<8} b={target:.4f}  sumZ={Z_t:7.2f}  rank {rank}/{len(bases)}")


if __name__ == "__main__":
    main()
