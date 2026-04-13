"""
Toroidal Proton: Electromagnetic Form Factor Test
==================================================

The definitive test of the "proton is a torus" hypothesis.

When electrons scatter off protons, the cross-section depends on the
electric form factor G_E(Q^2), which is the 3-D Fourier transform of
the charge distribution.  This has been measured with high precision.

We compute G_E(Q^2) for a toroidal charge distribution and compare it
directly to the experimental data.  No fitting of quark positions or
moments — just the shape of the torus vs what the scattering data says.

If no torus geometry matches the data, the hypothesis is ruled out.
If some geometry matches, those parameters are determined by the form
factor (not by moments), and the magnetic moment becomes a genuine
independent prediction.

The form factor of a charge distribution rho(r) is:

    G_E(Q^2) = (1/Q_total) * integral rho(r) * exp(i q.r) d^3r

For an axially symmetric distribution we can reduce this to a 2-D
integral using cylindrical coordinates and a Bessel function identity.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0   # zeroth-order Bessel function of the first kind

# =========================================================================
# EXPERIMENTAL DATA
# =========================================================================
# Proton electric form factor G_E(Q^2), normalised so G_E(0) = 1.
#
# The standard parameterisation is the "dipole form factor":
#     G_E_dipole(Q^2) = 1 / (1 + Q^2/0.71)^2
# where Q^2 is in GeV^2.  This fits the data well up to Q^2 ~ 5 GeV^2.
#
# More precise data from JLab polarisation transfer experiments show
# deviations at high Q^2, but the dipole is an excellent reference.

hbar_c = 0.1973  # GeV·fm  (conversion factor: hbar*c)

def G_E_dipole(Q2):
    """Standard dipole form factor (Q2 in GeV^2)."""
    return 1.0 / (1.0 + Q2 / 0.71)**2

# "Data points" from the dipole fit, which represents the experimental
# measurements well.  We also include the Kelly (2004) parameterisation
# which captures the high-Q^2 deviations seen at JLab.
def G_E_kelly(Q2):
    """Kelly (2004) parameterisation of proton G_E.
    Fits world data including JLab polarisation-transfer results."""
    tau = Q2 / (4 * 0.9383**2)   # tau = Q^2 / (4 * m_p^2)
    num = 1.0 - 0.24 * tau
    den = 1.0 + 10.98 * tau + 12.82 * tau**2 + 21.97 * tau**3
    return num / den

# =========================================================================
# TORUS CHARGE DISTRIBUTION AND ITS FORM FACTOR
# =========================================================================
#
# A torus with major radius R and minor radius r, centred at the origin,
# with symmetry axis along z.  We consider three charge distributions:
#
#   (a) Uniform volume charge:  rho = const inside the torus
#   (b) Uniform surface charge: rho on the torus surface only
#   (c) Gaussian profile:       rho falls off as exp(-d^2/2sigma^2)
#                                where d = distance from the torus skeleton
#
# The form factor for an axially symmetric distribution is:
#
#   G_E(q) = (1/Q_total) * 2*pi * integral rho(rho_cyl, z) *
#            J_0(q_perp * rho_cyl) * exp(i*q_z*z) * rho_cyl  d(rho_cyl) dz
#
# Since the distribution is symmetric under z -> -z reflection (for the
# torus), we can average over all orientations of q relative to the
# symmetry axis (the proton is spin-1/2 so in an unpolarised experiment
# we measure the orientation-averaged form factor).
#
# For a spin-1/2 particle, the form factor measured in elastic scattering
# is the SPHERICAL average:
#
#   G_E(Q^2) = <G_E(q)>_orientations
#            = (1/Q_tot) * integral rho(r) * sin(qr)/(qr) * d^3r
#
# where q = |q| = sqrt(Q^2 + ...) ~ sqrt(Q^2) in the Breit frame,
# and sin(qr)/(qr) is the angle-averaged plane-wave factor.
#
# We evaluate this numerically on a grid in cylindrical coordinates.

def torus_form_factor_volume(Q2_array, R, r, N_rho=200, N_z=200):
    """
    Compute the orientation-averaged electric form factor G_E(Q^2)
    for a torus with UNIFORM VOLUME charge.

    Parameters:
        Q2_array : array of Q^2 values in GeV^2
        R        : major radius in fm
        r        : minor radius in fm
        N_rho, N_z : grid resolution for numerical integration

    Returns:
        G_E : array of form factor values, normalised so G_E(0) = 1
    """
    # Set up cylindrical grid (rho_cyl, z)
    # The torus cross-section at cylindrical radius rho_cyl
    # exists where (rho_cyl - R)^2 + z^2 < r^2
    # So rho_cyl ranges from R-r to R+r, z from -r to r.

    rho_min = max(R - r, 1e-6)
    rho_max = R + r
    rho_arr = np.linspace(rho_min, rho_max, N_rho)
    z_arr = np.linspace(-r, r, N_z)

    d_rho = rho_arr[1] - rho_arr[0]
    d_z = z_arr[1] - z_arr[0]

    # Build the charge density on the grid (1 inside torus, 0 outside)
    RHO, Z = np.meshgrid(rho_arr, z_arr, indexing='ij')  # (N_rho, N_z)
    dist2 = (RHO - R)**2 + Z**2
    inside = dist2 <= r**2   # boolean mask

    # Radial distance from origin for each grid point: sqrt(rho_cyl^2 + z^2)
    r_3d = np.sqrt(RHO**2 + Z**2)

    # Integration weight: 2*pi*rho_cyl (from azimuthal integration) * d_rho * d_z
    # The factor 2*pi comes from integrating over the azimuthal angle phi.
    weight = 2 * np.pi * RHO * d_rho * d_z
    weight_inside = weight * inside   # zero outside the torus

    # Total charge (normalisation)
    Q_total = np.sum(weight_inside)

    # For each Q^2, compute G_E = (1/Q_total) * sum_grid [ weight * sin(qr)/(qr) ]
    G_E = np.zeros_like(Q2_array, dtype=float)
    for i, Q2 in enumerate(Q2_array):
        q = np.sqrt(Q2) / hbar_c   # momentum transfer in fm^-1
        qr = q * r_3d
        # sin(qr)/(qr), handling qr=0
        sinc = np.where(qr < 1e-12, 1.0, np.sin(qr) / qr)
        G_E[i] = np.sum(weight_inside * sinc) / Q_total

    return G_E


def torus_form_factor_shell(Q2_array, R, r, N_theta=500, N_phi=500):
    """
    Form factor for a torus with UNIFORM SURFACE charge.
    Parameterise the surface and integrate sin(qr)/(qr) over it.
    """
    theta = np.linspace(0, 2 * np.pi, N_theta, endpoint=False)  # toroidal angle
    phi   = np.linspace(0, 2 * np.pi, N_phi, endpoint=False)    # poloidal angle
    d_theta = theta[1] - theta[0]
    d_phi = phi[1] - phi[0]

    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')

    # Surface coordinates
    X = (R + r * np.cos(PHI)) * np.cos(THETA)
    Y = (R + r * np.cos(PHI)) * np.sin(THETA)
    Z = r * np.sin(PHI)

    # Distance from origin
    r_3d = np.sqrt(X**2 + Y**2 + Z**2)

    # Surface area element: |dr/dtheta x dr/dphi| d_theta d_phi
    # For a torus: dA = r * (R + r*cos(phi)) * d_theta * d_phi
    dA = r * (R + r * np.cos(PHI)) * d_theta * d_phi
    A_total = np.sum(dA)

    G_E = np.zeros_like(Q2_array, dtype=float)
    for i, Q2 in enumerate(Q2_array):
        q = np.sqrt(Q2) / hbar_c
        qr = q * r_3d
        sinc = np.where(qr < 1e-12, 1.0, np.sin(qr) / qr)
        G_E[i] = np.sum(dA * sinc) / A_total

    return G_E


def torus_form_factor_gaussian(Q2_array, R, sigma, N_rho=200, N_z=200):
    """
    Form factor for a torus with GAUSSIAN charge profile.
    Charge density ~ exp(-d^2 / (2*sigma^2)) where d is the distance
    from the torus skeleton (the central circle of radius R).
    sigma controls the tube thickness (effective minor radius ~ 2*sigma).
    """
    r_eff = 4 * sigma   # integrate out to 4 sigma
    rho_min = max(R - r_eff, 1e-6)
    rho_max = R + r_eff
    rho_arr = np.linspace(rho_min, rho_max, N_rho)
    z_arr = np.linspace(-r_eff, r_eff, N_z)

    d_rho = rho_arr[1] - rho_arr[0]
    d_z = z_arr[1] - z_arr[0]

    RHO, Z = np.meshgrid(rho_arr, z_arr, indexing='ij')
    dist2 = (RHO - R)**2 + Z**2   # distance^2 from torus skeleton
    density = np.exp(-dist2 / (2 * sigma**2))

    r_3d = np.sqrt(RHO**2 + Z**2)
    weight = 2 * np.pi * RHO * d_rho * d_z * density
    Q_total = np.sum(weight)

    G_E = np.zeros_like(Q2_array, dtype=float)
    for i, Q2 in enumerate(Q2_array):
        q = np.sqrt(Q2) / hbar_c
        qr = q * r_3d
        sinc = np.where(qr < 1e-12, 1.0, np.sin(qr) / qr)
        G_E[i] = np.sum(weight * sinc) / Q_total

    return G_E


# =========================================================================
# ALSO COMPUTE THE FORM FACTOR FOR A GAUSSIAN SPHERE (for comparison)
# =========================================================================

def sphere_form_factor_gaussian(Q2_array, r_rms):
    """
    Form factor for a Gaussian charge distribution (spherical):
        rho(r) ~ exp(-r^2 / (2*sigma^2)),  sigma = r_rms / sqrt(3)
    This has the exact analytical form factor:
        G_E = exp(-Q^2 * r_rms^2 / (6 * hbar_c^2))
    """
    return np.exp(-Q2_array * r_rms**2 / (6 * hbar_c**2))


# =========================================================================
# COMPUTE AND COMPARE
# =========================================================================

Q2 = np.linspace(0, 5.0, 300)   # Q^2 in GeV^2

# Experimental references
G_dip = G_E_dipole(Q2)
G_kel = G_E_kelly(Q2)

# Gaussian sphere (for reference)
G_sphere = sphere_form_factor_gaussian(Q2, 0.841)

print("=" * 70)
print("  TOROIDAL PROTON: FORM FACTOR TEST")
print("=" * 70)
print("\nComputing form factors (this may take a moment)...")

# --- Scan over torus geometries ---
# We try several (R, r) combinations.
# The charge radius of a uniform-volume torus is:
#     <r^2> = R^2 + 3r^2/4 + r^4/(4R^2)  [approximately R^2 + 3r^2/4]
# We choose combinations that give roughly the right charge radius (~0.84 fm).

torus_configs = [
    # (R, r, label, description)
    (0.70, 0.30, "fat",    "R=0.70, r=0.30 fm (R/r=2.3, fat torus)"),
    (0.60, 0.45, "donut",  "R=0.60, r=0.45 fm (R/r=1.3, nearly sphere-like)"),
    (0.75, 0.25, "mod",    "R=0.75, r=0.25 fm (R/r=3.0, moderate)"),
    (0.80, 0.15, "thin",   "R=0.80, r=0.15 fm (R/r=5.3, thin ring)"),
    (0.50, 0.50, "horn",   "R=0.50, r=0.50 fm (R/r=1.0, horn torus)"),
]

# Compute charge radii for each configuration
print(f"\n  {'Config':<8} {'R [fm]':>8} {'r [fm]':>8} {'R/r':>6} {'r_ch [fm]':>10}")
print(f"  {'-'*44}")

torus_results = {}
for R, r, label, desc in torus_configs:
    # Approximate RMS charge radius for uniform-volume torus
    r_ch = np.sqrt(R**2 + 3*r**2/4)
    print(f"  {label:<8} {R:>8.2f} {r:>8.2f} {R/r:>6.1f} {r_ch:>10.3f}")

    # Compute form factor (uniform volume)
    G_vol = torus_form_factor_volume(Q2, R, r, N_rho=250, N_z=250)
    torus_results[label] = {
        'R': R, 'r': r, 'G_vol': G_vol,
        'desc': desc, 'r_ch': r_ch
    }

# Also compute Gaussian torus with best-guess parameters
print("\nComputing Gaussian torus profiles...")
gauss_configs = [
    (0.65, 0.25, "gauss1", "Gaussian: R=0.65, sigma=0.25 fm"),
    (0.55, 0.35, "gauss2", "Gaussian: R=0.55, sigma=0.35 fm"),
]
for R, sigma, label, desc in gauss_configs:
    r_ch = np.sqrt(R**2 + 3*sigma**2)
    print(f"  {label:<8}  R={R:.2f}, sigma={sigma:.2f}  ->  r_ch = {r_ch:.3f} fm")
    G_gauss = torus_form_factor_gaussian(Q2, R, sigma, N_rho=250, N_z=250)
    torus_results[label] = {
        'R': R, 'r': sigma, 'G_vol': G_gauss,
        'desc': desc, 'r_ch': r_ch
    }


# =========================================================================
# FIT QUALITY: find which torus best matches the data
# =========================================================================
print(f"\n  Fit quality (chi^2-like, lower = better match to dipole):")
print(f"  {'Config':<8} {'sum |G_torus - G_dipole|^2':>30}")
print(f"  {'-'*40}")

best_label = None
best_chi2 = 1e10
for label, res in torus_results.items():
    chi2 = np.sum((res['G_vol'] - G_dip)**2) * (Q2[1] - Q2[0])
    res['chi2'] = chi2
    if chi2 < best_chi2:
        best_chi2 = chi2
        best_label = label
    print(f"  {label:<8} {chi2:>30.6f}")

# Gaussian sphere chi2
chi2_sphere = np.sum((G_sphere - G_dip)**2) * (Q2[1] - Q2[0])
print(f"  {'sphere':<8} {chi2_sphere:>30.6f}  (Gaussian sphere, r_rms=0.841)")

print(f"\n  Best torus: {best_label} ({torus_results[best_label]['desc']})")


# =========================================================================
# KEY DIAGNOSTIC: does the torus form factor have the right SHAPE?
# =========================================================================
# The dipole falls off as 1/Q^4 at large Q.  A torus has oscillations
# (diffraction minima) because it has a hole in the middle.  A sphere
# doesn't.  These oscillations would be the smoking gun.

print(f"\n{'='*70}")
print("  KEY QUESTION: Does a torus produce diffraction minima?")
print(f"{'='*70}")

# Check if any torus form factor crosses zero
for label, res in torus_results.items():
    G = res['G_vol']
    zero_crossings = np.where(np.diff(np.sign(G)))[0]
    if len(zero_crossings) > 0:
        Q2_zeros = Q2[zero_crossings]
        print(f"  {label}: G_E crosses zero at Q^2 = {', '.join(f'{z:.2f}' for z in Q2_zeros[:3])} GeV^2")
    else:
        print(f"  {label}: G_E does NOT cross zero in 0-5 GeV^2")

print(f"\n  Experimental G_E: no confirmed zero crossing below Q^2 ~ 5 GeV^2")
print(f"  (Some analyses suggest a possible zero near Q^2 ~ 8-10 GeV^2)")


# =========================================================================
# COMPUTE MAGNETIC MOMENT FROM BEST-FIT TORUS (the non-circular test)
# =========================================================================

print(f"\n{'='*70}")
print("  NON-CIRCULAR PREDICTION: magnetic moment from form-factor geometry")
print(f"{'='*70}")

# For the best-fitting torus, compute what the magnetic moment would be
# if we use the form-factor-determined geometry.
best = torus_results[best_label]
R_best, r_best = best['R'], best['r']

# Simple circular orbit model: mu = e*v*R / 2
# In nuclear magnetons: mu/mu_N = m_p * c * R / hbar = R / (hbar/(m_p*c))
# = R / 0.2103 fm (proton Compton wavelength)
lambda_p = 0.2103   # proton Compton wavelength in fm

# If the charge simply orbits at radius R (like a spinning top):
mu_simple = R_best / lambda_p   # in nuclear magnetons
print(f"\n  Best-fit torus: R = {R_best:.2f} fm, r = {r_best:.2f} fm")
print(f"  If charge orbits at R:  mu = {mu_simple:.2f} mu_N  (experiment: 2.79)")

# More carefully: for a uniform volume torus spinning about its axis,
# mu = Q * omega * <rho^2> / 2, L = M * omega * <rho^2>
# mu/L = Q/(2M), so mu = (Q/2M) * L = mu_N * (L/hbar)
# For spin-1/2: L = hbar/2, so mu = mu_N/2 = 0.5 mu_N (Dirac value)
# The anomalous part requires non-uniform charge/mass distribution.
print(f"  Uniform spinning torus:  mu = 0.50 mu_N  (Dirac value for any shape)")
print(f"  The anomalous moment (2.79 vs 0.50) requires charge != mass distribution")
print(f"  This is where internal structure (quarks/flux tubes) re-enters")


# =========================================================================
# VISUALIZATION
# =========================================================================

fig, axes = plt.subplots(1, 3, figsize=(20, 6))

# --- Panel 1: Form factors on linear scale ---
ax1 = axes[0]
ax1.plot(Q2, G_dip, 'k-', linewidth=2.5, label='Experiment (dipole)')
ax1.plot(Q2, G_kel, 'k--', linewidth=1.5, label='Experiment (Kelly fit)')
ax1.plot(Q2, G_sphere, 'gray', linewidth=1.5, linestyle=':', label='Gaussian sphere')

colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628', '#f781bf']
for i, (label, res) in enumerate(torus_results.items()):
    ax1.plot(Q2, res['G_vol'], color=colors[i % len(colors)],
             linewidth=1.5, alpha=0.8, label=res['desc'][:30])

ax1.set_xlabel('Q² [GeV²]', fontsize=12)
ax1.set_ylabel('G_E(Q²)', fontsize=12)
ax1.set_title('Form Factor (linear scale)', fontsize=12)
ax1.legend(fontsize=7, loc='upper right')
ax1.set_ylim(-0.15, 1.05)
ax1.axhline(y=0, color='gray', linewidth=0.5)
ax1.grid(True, alpha=0.3)

# --- Panel 2: Form factor ratio to dipole ---
ax2 = axes[1]
ax2.axhline(y=1, color='black', linewidth=2, label='Dipole (reference)')

mask = Q2 > 0.1   # avoid division by zero region
for i, (label, res) in enumerate(torus_results.items()):
    ratio = res['G_vol'][mask] / G_dip[mask]
    ax2.plot(Q2[mask], ratio, color=colors[i % len(colors)],
             linewidth=1.5, alpha=0.8, label=label)

ratio_kelly = G_kel[mask] / G_dip[mask]
ax2.plot(Q2[mask], ratio_kelly, 'k--', linewidth=1.5, label='Kelly/dipole')

ax2.set_xlabel('Q² [GeV²]', fontsize=12)
ax2.set_ylabel('G_E / G_E(dipole)', fontsize=12)
ax2.set_title('Ratio to dipole', fontsize=12)
ax2.legend(fontsize=8)
ax2.set_ylim(0, 2.5)
ax2.grid(True, alpha=0.3)

# --- Panel 3: Log scale to see high-Q^2 behaviour ---
ax3 = axes[2]
ax3.semilogy(Q2[1:], np.abs(G_dip[1:]), 'k-', linewidth=2.5, label='Experiment (dipole)')
ax3.semilogy(Q2[1:], np.abs(G_kel[1:]), 'k--', linewidth=1.5, label='Kelly fit')
ax3.semilogy(Q2[1:], np.abs(G_sphere[1:]), 'gray', linewidth=1.5, linestyle=':',
             label='Gaussian sphere')

for i, (label, res) in enumerate(torus_results.items()):
    G = res['G_vol']
    ax3.semilogy(Q2[1:], np.abs(G[1:]), color=colors[i % len(colors)],
                 linewidth=1.5, alpha=0.8, label=label)

ax3.set_xlabel('Q² [GeV²]', fontsize=12)
ax3.set_ylabel('|G_E(Q²)|', fontsize=12)
ax3.set_title('Form factor (log scale)', fontsize=12)
ax3.legend(fontsize=8)
ax3.set_ylim(1e-4, 1.5)
ax3.grid(True, alpha=0.3)

fig.suptitle('Proton Form Factor: Torus vs Experiment', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('toroidal_proton_form_factor.png', dpi=150, bbox_inches='tight')
print(f"\nPlot saved to toroidal_proton_form_factor.png")
plt.show()
