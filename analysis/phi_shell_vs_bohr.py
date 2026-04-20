"""
Direct test of 'shell radii follow r_n = r_0 * phi^n'.

Bohr/Schrödinger (observed): E_n proportional to -1/n^2.
Phi-shell hypothesis:        E_n proportional to -1/phi^n.

For Lyman series (transitions n -> 1), the photon energy equals |E_1| - |E_n|:
  observed:      E_gamma(n) = R * (1 - 1/n^2)
  phi-shell:     E_gamma(n) = C * (1 - 1/phi^(n-1))   [normalized so E_gamma(inf)=C]

Compare predicted wavelengths to NIST values for Lyman-alpha..epsilon.
"""

import numpy as np

PHI = (1 + 5 ** 0.5) / 2
R_INF_NM = 91.1267  # Rydberg series limit (vacuum), nm (hc/R ~ 91.127 nm)

# NIST Lyman-series vacuum wavelengths (nm), n = 2..6
LYMAN = {
    2: 121.567,
    3: 102.572,
    4:  97.253,
    5:  94.974,
    6:  93.780,
}

# (1) Bohr 1/n^2 model:
#     lambda_n = R_inf_nm / (1 - 1/n^2)
# (2) Phi-shell 1/phi^n model:
#     lambda_n = (R_inf_nm / some factor) / (1 - 1/phi^(n-1))
# Normalize each model to hit Ly-alpha exactly (one-parameter anchor, just
# like the cosmology paper anchors its comoving distance to the CMB).

def fit_anchor(n_anchor, lam_obs, model_shape):
    """model_shape(n) gives the dimensionless 1 - 1/f(n) factor.
    Return the constant C such that lam_pred(n) = C / model_shape(n)
    matches lam_obs at n=n_anchor."""
    return lam_obs * model_shape(n_anchor)


def bohr_shape(n):
    return 1.0 - 1.0 / n ** 2


def phi_shape(n):
    # Shells indexed from n=1 ground state. Photon energy ~ 1 - 1/phi^(n-1).
    return 1.0 - 1.0 / PHI ** (n - 1)


def main():
    n_vals = sorted(LYMAN.keys())
    lam_obs = np.array([LYMAN[n] for n in n_vals])

    C_bohr = fit_anchor(n_vals[0], lam_obs[0], bohr_shape)
    C_phi  = fit_anchor(n_vals[0], lam_obs[0], phi_shape)

    print(f"Anchored both models to Ly-alpha (n=2).")
    print(f"{'n':>3} {'observed':>12} {'Bohr pred':>12} {'err %':>9}"
          f"  {'phi pred':>12} {'err %':>9}")
    for n, lo in zip(n_vals, lam_obs):
        lb = C_bohr / bohr_shape(n)
        lp = C_phi / phi_shape(n)
        print(f"{n:>3} {lo:>12.4f} {lb:>12.4f} {(lb/lo-1)*100:>+9.3f}"
              f"  {lp:>12.4f} {(lp/lo-1)*100:>+9.3f}")

    # RMS relative error, excluding the anchor point
    rms_bohr = np.sqrt(np.mean([((C_bohr/bohr_shape(n))/LYMAN[n]-1)**2
                                for n in n_vals[1:]])) * 100
    rms_phi  = np.sqrt(np.mean([((C_phi/phi_shape(n))/LYMAN[n]-1)**2
                                for n in n_vals[1:]])) * 100
    print(f"\nRMS relative error (excluding anchor):")
    print(f"  Bohr 1/n^2 model:     {rms_bohr:.4f} %")
    print(f"  phi-shell 1/phi^n:    {rms_phi:.4f} %")
    print(f"  ratio:                {rms_phi/rms_bohr:.1f}x worse" if rms_bohr > 0
          else "")


if __name__ == "__main__":
    main()
