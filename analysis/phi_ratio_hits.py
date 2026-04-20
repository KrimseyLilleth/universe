"""
Pure ratio-inspection: for each element, look at all pairwise wavelength
and frequency ratios and report which (if any) fall close to phi^k for
k = -4..4.

No model is imposed. This is just counting hits in existing published
spectral data.
"""

import numpy as np
from spectral_phi_analysis import ELEMENTS, pairwise_ratios

PHI = (1 + 5 ** 0.5) / 2


def closest_phi_power(r, k_range=range(-4, 5)):
    """Return (best_k, fractional_error) such that r is closest to phi^best_k."""
    best_k, best_err = None, float("inf")
    for k in k_range:
        target = PHI ** k
        err = abs(r / target - 1.0)
        if err < best_err:
            best_err = err
            best_k = k
    return best_k, best_err


def report_element(name, lines, tol=0.005):
    """List pairs whose frequency ratio is within `tol` of some phi^k (k!=0)."""
    labels = list(lines.keys())
    wls = np.array([lines[k] for k in labels])
    nu = 1.0 / wls

    print(f"\n=== {name} ({len(labels)} lines, {len(labels)*(len(labels)-1)//2} pairs) ===")
    print(f"Pairs with frequency ratio within {tol*100:.1f}% of phi^k (k in -4..4, k!=0):")
    hits = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            r = nu[i] / nu[j]
            if r < 1.0:
                r = 1.0 / r
                lo, hi = labels[i], labels[j]  # hi has larger wl
            else:
                lo, hi = labels[j], labels[i]
            k, err = closest_phi_power(r)
            if k != 0 and err < tol:
                hits.append((err, k, hi, lo, r))
    hits.sort()
    if not hits:
        print("  (none)")
    else:
        for err, k, hi, lo, r in hits:
            print(f"  {hi:<14} / {lo:<14}  ratio = {r:7.4f}"
                  f"  phi^{k:+d} = {PHI**k:7.4f}  err = {err*100:+.3f}%")
    return hits


def expected_hits_empirical(ratios, k_range, tol, n_shuffles=500, rng=None):
    """Empirical null: randomly perturb each log-ratio by a uniform offset
    within a few log(phi) widths, count hits, average. This preserves the
    density of log-ratios but removes any quantum-mechanical structure."""
    if rng is None:
        rng = np.random.default_rng(0)
    logs = np.log(np.asarray(ratios))
    log_phi = np.log(PHI)
    targets_log = np.array([k * log_phi for k in k_range if k != 0])
    tol_log = np.log(1 + tol)
    counts = []
    for _ in range(n_shuffles):
        shift = rng.uniform(-log_phi / 2, log_phi / 2, size=logs.shape)
        shifted = logs + shift
        hit = np.zeros_like(shifted, dtype=bool)
        for t in targets_log:
            hit |= np.abs(shifted - t) < tol_log
        counts.append(int(hit.sum()))
    return float(np.mean(counts)), float(np.std(counts))


def main():
    print(f"phi   = {PHI:.6f}")
    print(f"phi^2 = {PHI**2:.6f}    phi^-1 = {1/PHI:.6f}")
    print(f"phi^3 = {PHI**3:.6f}    phi^-2 = {1/PHI**2:.6f}")
    print(f"phi^4 = {PHI**4:.6f}    phi^-3 = {1/PHI**3:.6f}")

    k_range = range(-4, 5)
    for tol in [0.005, 0.01, 0.02]:
        print(f"\n\n#################  tolerance = {tol*100:.1f}%  #################")
        total_hits = 0
        total_pairs = 0
        all_ratios = []
        for name, lines in ELEMENTS.items():
            hits = report_element(name, lines, tol=tol)
            total_hits += len(hits)
            total_pairs += len(lines) * (len(lines) - 1) // 2
            all_ratios.extend(pairwise_ratios(list(lines.values())).tolist())
        mean_expected, std_expected = expected_hits_empirical(
            all_ratios, k_range, tol)
        z = (total_hits - mean_expected) / std_expected if std_expected > 0 else 0
        print(f"\nTotal hits: {total_hits} out of {total_pairs} pairs")
        print(f"Empirical null (random log-shift of same ratios, 500 trials):")
        print(f"  expected = {mean_expected:.1f} +/- {std_expected:.1f}")
        print(f"  observed z = {z:+.2f}")


if __name__ == "__main__":
    main()
