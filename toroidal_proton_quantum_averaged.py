"""
Quantum-Averaged Toroidal Proton: Can fluctuations hide the hole?
=================================================================

The rigid torus form factor has diffraction zeros that experiment
doesn't see.  But a quantum proton isn't rigid — its internal geometry
fluctuates.  If the torus size fluctuates quantum mechanically, different
configurations have zeros at different Q^2, and the average might wash
them out.

We test this with two approaches:
  (A) Quantum size fluctuations: average the torus form factor over a
      Gaussian distribution of major radii R.
  (B) Torus core + spherical halo: a fraction of the charge is toroidal,
      the rest is in a spherical cloud (gluon/pion cloud).

The question: can either approach reproduce the measured form factor?
And if so, how much toroidal structure actually survives?
"""

import numpy as np
import matplotlib.pyplot as plt

# =========================================================================
# SETUP
# =========================================================================

hbar_c = 0.1973    # GeV·fm

def G_E_dipole(Q2):
    """Standard dipole form factor (experimental reference)."""
    return 1.0 / (1.0 + Q2 / 0.71)**2

def G_E_kelly(Q2):
    """Kelly (2004) parameterisation of proton G_E."""
    tau = Q2 / (4 * 0.9383**2)
    return (1.0 - 0.24 * tau) / (1.0 + 10.98 * tau + 12.82 * tau**2 + 21.97 * tau**3)

# Q^2 grid
N_Q = 200
Q2 = np.linspace(0, 5.0, N_Q)

# Spatial grid in cylindrical coordinates (rho_cyl, z)
# This is the cross-section plane — we exploit azimuthal symmetry.
N_grid = 200
rho_arr = np.linspace(0.01, 3.0, N_grid)   # cylindrical radius [fm]
z_arr   = np.linspace(-2.5, 2.5, N_grid)   # height [fm]
d_rho = rho_arr[1] - rho_arr[0]
d_z   = z_arr[1]   - z_arr[0]

RHO, Z = np.meshgrid(rho_arr, z_arr, indexing='ij')   # (N_grid, N_grid)
r_3d   = np.sqrt(RHO**2 + Z**2)       # distance from origin

# Azimuthal integration weight: 2*pi*rho * d_rho * d_z
az_weight = 2.0 * np.pi * RHO * d_rho * d_z

# =========================================================================
# PRECOMPUTE sin(q*r)/(q*r) FOR EVERY (Q^2, grid point) — the expensive part
# =========================================================================
print("Precomputing sinc grid...")

q_arr = np.sqrt(Q2) / hbar_c     # momentum transfer in fm^-1
# Shape: (N_Q, N_grid, N_grid)
sinc_grid = np.zeros((N_Q, N_grid, N_grid))
for i, q in enumerate(q_arr):
    qr = q * r_3d
    sinc_grid[i] = np.where(qr < 1e-12, 1.0, np.sin(qr) / qr)

print("  Done.")


def gaussian_torus_density(R, sigma):
    """Charge density for a Gaussian torus: exp(-d^2 / (2*sigma^2))
    where d = distance from the circle of radius R in the z=0 plane."""
    dist2 = (RHO - R)**2 + Z**2
    return np.exp(-dist2 / (2.0 * sigma**2))

def gaussian_sphere_density(sigma_s):
    """Spherical Gaussian charge density: exp(-r^2 / (2*sigma_s^2))."""
    return np.exp(-r_3d**2 / (2.0 * sigma_s**2))

def form_factor_from_density(density):
    """Compute G_E(Q^2) from a density on the precomputed grid.
    Uses the precomputed sinc_grid for speed."""
    weighted = density * az_weight
    Q_total = np.sum(weighted)
    if Q_total < 1e-30:
        return np.ones(N_Q)
    # G_E[i] = sum(weighted * sinc_grid[i]) / Q_total
    return np.einsum('jk,ijk->i', weighted, sinc_grid) / Q_total


# =========================================================================
# APPROACH A: QUANTUM SIZE FLUCTUATIONS
# =========================================================================
# Average the Gaussian-torus form factor over a distribution of R values.
# The tube thickness sigma is held fixed.
# Different R values put the diffraction zeros at different Q^2,
# so the average may fill them in.

print("\n" + "=" * 70)
print("  APPROACH A: QUANTUM SIZE FLUCTUATIONS")
print("=" * 70)

# Torus parameters
R0 = 0.65          # mean major radius [fm]
sigma_tube = 0.30  # tube thickness (Gaussian width) [fm]

# Fluctuation widths to try
sigma_R_values = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]

print(f"\n  Mean torus: R0 = {R0} fm, tube sigma = {sigma_tube} fm")
print(f"\n  {'sigma_R [fm]':>12}  {'sigma_R/R0':>10}  {'has zeros?':>12}  {'chi2 vs dipole':>16}")
print(f"  {'-'*55}")

results_A = {}
for sigma_R in sigma_R_values:
    if sigma_R < 0.001:
        # No fluctuation: just the single torus
        rho_density = gaussian_torus_density(R0, sigma_tube)
        G_avg = form_factor_from_density(rho_density)
    else:
        # Sample R from Gaussian and average
        N_samples = 50
        R_lo = max(0.01, R0 - 3.5 * sigma_R)
        R_hi = R0 + 3.5 * sigma_R
        R_samples = np.linspace(R_lo, R_hi, N_samples)
        dR = R_samples[1] - R_samples[0]

        G_avg = np.zeros(N_Q)
        total_w = 0
        for R_s in R_samples:
            w = np.exp(-(R_s - R0)**2 / (2 * sigma_R**2)) * dR
            rho_s = gaussian_torus_density(R_s, sigma_tube)
            G_avg += w * form_factor_from_density(rho_s)
            total_w += w
        G_avg /= total_w

    # Check for zeros
    zero_crossings = np.where(np.diff(np.sign(G_avg)))[0]
    has_zeros = len(zero_crossings) > 0

    # Fit quality
    chi2 = np.sum((G_avg - G_E_dipole(Q2))**2) * (Q2[1] - Q2[0])

    frac = sigma_R / R0 if R0 > 0 else 0
    zeros_str = "YES" if has_zeros else "NO"
    print(f"  {sigma_R:>12.2f}  {frac:>10.1%}  {zeros_str:>12}  {chi2:>16.6f}")

    results_A[sigma_R] = {'G': G_avg, 'chi2': chi2, 'has_zeros': has_zeros}


# =========================================================================
# APPROACH B: TORUS CORE + SPHERICAL HALO
# =========================================================================
# Model: rho = alpha * rho_torus + (1-alpha) * rho_sphere
# The torus is the "valence" structure; the sphere is the pion/gluon cloud.

print(f"\n{'='*70}")
print("  APPROACH B: TORUS CORE + SPHERICAL HALO")
print(f"{'='*70}")

# Torus core
R_core = 0.55
sigma_core = 0.25
rho_torus = gaussian_torus_density(R_core, sigma_core)
G_torus = form_factor_from_density(rho_torus)

# Spherical halo (RMS radius ~ 0.84 fm -> sigma = 0.84/sqrt(3) ~ 0.49 fm)
sigma_halo = 0.50
rho_sphere = gaussian_sphere_density(sigma_halo)
G_sphere = form_factor_from_density(rho_sphere)

alpha_values = np.arange(0.0, 1.05, 0.05)

print(f"\n  Torus core: R = {R_core}, sigma = {sigma_core} fm")
print(f"  Spherical halo: sigma = {sigma_halo} fm")
print(f"\n  {'alpha':>8}  {'has zeros?':>12}  {'chi2 vs dipole':>16}")
print(f"  {'-'*40}")

results_B = {}
best_alpha = 0
best_chi2_B = 1e10
for alpha in alpha_values:
    G_mix = alpha * G_torus + (1 - alpha) * G_sphere
    zero_crossings = np.where(np.diff(np.sign(G_mix)))[0]
    has_zeros = len(zero_crossings) > 0
    chi2 = np.sum((G_mix - G_E_dipole(Q2))**2) * (Q2[1] - Q2[0])

    if chi2 < best_chi2_B:
        best_chi2_B = chi2
        best_alpha = alpha

    results_B[alpha] = {'G': G_mix, 'chi2': chi2, 'has_zeros': has_zeros}
    if alpha * 100 % 10 < 1:  # print every 10%
        zs = "YES" if has_zeros else "NO"
        print(f"  {alpha:>8.2f}  {zs:>12}  {chi2:>16.6f}")

# Find the maximum alpha that avoids zeros
max_alpha_no_zeros = 0
for alpha in sorted(results_B.keys()):
    if not results_B[alpha]['has_zeros']:
        max_alpha_no_zeros = alpha

print(f"\n  Best fit: alpha = {best_alpha:.2f}  (chi2 = {best_chi2_B:.6f})")
print(f"  Maximum toroidal fraction without zeros: alpha = {max_alpha_no_zeros:.2f}")
print(f"  Meaning: at most {max_alpha_no_zeros*100:.0f}% of the charge can be toroidal")


# =========================================================================
# APPROACH C: COMBINED — size fluctuations + halo
# =========================================================================

print(f"\n{'='*70}")
print("  APPROACH C: FLUCTUATING TORUS + SPHERICAL HALO (combined)")
print(f"{'='*70}")

# Use the best fluctuation width from approach A and mix with sphere
sigma_R_use = 0.30  # significant quantum fluctuation
N_samp = 50
R_lo = max(0.01, R0 - 3.5 * sigma_R_use)
R_hi = R0 + 3.5 * sigma_R_use
R_samps = np.linspace(R_lo, R_hi, N_samp)
dR = R_samps[1] - R_samps[0]

G_fluct_torus = np.zeros(N_Q)
total_w = 0
for Rs in R_samps:
    w = np.exp(-(Rs - R0)**2 / (2 * sigma_R_use**2)) * dR
    rho_s = gaussian_torus_density(Rs, sigma_tube)
    G_fluct_torus += w * form_factor_from_density(rho_s)
    total_w += w
G_fluct_torus /= total_w

print(f"\n  Fluctuating torus: R0={R0}, sigma_R={sigma_R_use}, tube_sigma={sigma_tube}")
print(f"  Spherical halo: sigma={sigma_halo}")
print(f"\n  {'alpha':>8}  {'has zeros?':>12}  {'chi2':>16}")
print(f"  {'-'*40}")

best_chi2_C = 1e10
best_alpha_C = 0
results_C = {}
for alpha in alpha_values:
    G_mix = alpha * G_fluct_torus + (1 - alpha) * G_sphere
    zero_crossings = np.where(np.diff(np.sign(G_mix)))[0]
    has_zeros = len(zero_crossings) > 0
    chi2 = np.sum((G_mix - G_E_dipole(Q2))**2) * (Q2[1] - Q2[0])
    results_C[alpha] = {'G': G_mix, 'chi2': chi2, 'has_zeros': has_zeros}
    if chi2 < best_chi2_C:
        best_chi2_C = chi2
        best_alpha_C = alpha
    if alpha * 100 % 10 < 1:
        zs = "YES" if has_zeros else "NO"
        print(f"  {alpha:>8.2f}  {zs:>12}  {chi2:>16.6f}")

max_alpha_C = 0
for alpha in sorted(results_C.keys()):
    if not results_C[alpha]['has_zeros']:
        max_alpha_C = alpha

print(f"\n  Best fit: alpha = {best_alpha_C:.2f}  (chi2 = {best_chi2_C:.6f})")
print(f"  Maximum toroidal fraction without zeros: alpha = {max_alpha_C:.2f}")


# =========================================================================
# SUMMARY
# =========================================================================

print(f"\n{'='*70}")
print("  SUMMARY")
print(f"{'='*70}")

print(f"""
  The question: can quantum averaging hide a toroidal proton?

  Approach A (size fluctuations only):
    sigma_R must be >= {min(s for s, r in results_A.items() if not r['has_zeros']) if any(not r['has_zeros'] for r in results_A.values()) else 'N/A'} fm
    to eliminate diffraction zeros.
    That's {min(s for s, r in results_A.items() if not r['has_zeros'])/R0*100 if any(not r['has_zeros'] for r in results_A.values()) else 0:.0f}% fluctuation in R — {'modest' if min((s for s, r in results_A.items() if not r['has_zeros']), default=999)/R0 < 0.3 else 'very large'}.

  Approach B (torus + sphere mix, no fluctuations):
    At most {max_alpha_no_zeros*100:.0f}% of the charge can be toroidal.
    Best fit at alpha = {best_alpha:.0%} toroidal.

  Approach C (fluctuating torus + sphere):
    Best fit at alpha = {best_alpha_C:.0%} toroidal.
    Maximum toroidal fraction without zeros: {max_alpha_C:.0%}.
""")


# =========================================================================
# VISUALIZATION
# =========================================================================

fig = plt.figure(figsize=(20, 12))

# --- Panel 1: Approach A — form factors for different sigma_R ---
ax1 = fig.add_subplot(231)
ax1.plot(Q2, G_E_dipole(Q2), 'k-', linewidth=2.5, label='Experiment')
for sigma_R in [0.0, 0.10, 0.20, 0.30, 0.50]:
    if sigma_R in results_A:
        lbl = f'sigma_R = {sigma_R:.2f} fm'
        ax1.plot(Q2, results_A[sigma_R]['G'], linewidth=1.5, alpha=0.8, label=lbl)
ax1.axhline(0, color='gray', linewidth=0.5)
ax1.set_xlabel('Q² [GeV²]')
ax1.set_ylabel('G_E(Q²)')
ax1.set_title('(A) Size fluctuations')
ax1.legend(fontsize=7)
ax1.set_ylim(-0.15, 1.05)
ax1.grid(True, alpha=0.3)

# --- Panel 2: Approach B — form factors for different alpha ---
ax2 = fig.add_subplot(232)
ax2.plot(Q2, G_E_dipole(Q2), 'k-', linewidth=2.5, label='Experiment')
for alpha in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
    if alpha in results_B:
        lbl = f'alpha = {alpha:.1f}'
        ax2.plot(Q2, results_B[alpha]['G'], linewidth=1.5, alpha=0.8, label=lbl)
ax2.axhline(0, color='gray', linewidth=0.5)
ax2.set_xlabel('Q² [GeV²]')
ax2.set_ylabel('G_E(Q²)')
ax2.set_title(f'(B) Torus + sphere mix')
ax2.legend(fontsize=7)
ax2.set_ylim(-0.15, 1.05)
ax2.grid(True, alpha=0.3)

# --- Panel 3: Best fits compared ---
ax3 = fig.add_subplot(233)
ax3.plot(Q2, G_E_dipole(Q2), 'k-', linewidth=2.5, label='Experiment (dipole)')
ax3.plot(Q2, G_E_kelly(Q2), 'k--', linewidth=1.5, label='Kelly fit')

# Best from each approach
best_sigR = min(results_A.keys(), key=lambda s: results_A[s]['chi2'])
ax3.plot(Q2, results_A[best_sigR]['G'], 'r-', linewidth=2,
         label=f'Best A: sigma_R={best_sigR:.2f}')
ax3.plot(Q2, results_B[best_alpha]['G'], 'b-', linewidth=2,
         label=f'Best B: alpha={best_alpha:.2f}')
ax3.plot(Q2, results_C[best_alpha_C]['G'], 'g-', linewidth=2,
         label=f'Best C: alpha={best_alpha_C:.2f}')

ax3.axhline(0, color='gray', linewidth=0.5)
ax3.set_xlabel('Q² [GeV²]')
ax3.set_ylabel('G_E(Q²)')
ax3.set_title('Best fits compared')
ax3.legend(fontsize=7)
ax3.set_ylim(-0.15, 1.05)
ax3.grid(True, alpha=0.3)

# --- Panel 4: Charge density cross-section for sharp torus ---
ax4 = fig.add_subplot(234)
rho_sharp = gaussian_torus_density(R0, sigma_tube)
im4 = ax4.pcolormesh(rho_arr, z_arr, rho_sharp.T, cmap='inferno', shading='auto')
ax4.set_xlabel('rho [fm]')
ax4.set_ylabel('z [fm]')
ax4.set_title(f'Sharp torus (R={R0}, sigma={sigma_tube})')
ax4.set_aspect('equal')
plt.colorbar(im4, ax=ax4, shrink=0.7, label='charge density')

# --- Panel 5: Quantum-averaged charge density ---
ax5 = fig.add_subplot(235)
# Build the averaged density
rho_avg = np.zeros_like(rho_sharp)
total_w = 0
sigma_R_best = best_sigR
for Rs in R_samps:
    w = np.exp(-(Rs - R0)**2 / (2 * sigma_R_best**2)) * dR
    rho_avg += w * gaussian_torus_density(Rs, sigma_tube)
    total_w += w
rho_avg /= total_w

im5 = ax5.pcolormesh(rho_arr, z_arr, rho_avg.T, cmap='inferno', shading='auto')
ax5.set_xlabel('rho [fm]')
ax5.set_ylabel('z [fm]')
ax5.set_title(f'Fluctuating torus (sigma_R={sigma_R_best:.2f})')
ax5.set_aspect('equal')
plt.colorbar(im5, ax=ax5, shrink=0.7, label='charge density')

# --- Panel 6: Mixture density (best alpha from C) ---
ax6 = fig.add_subplot(236)
rho_mix = best_alpha_C * rho_avg + (1 - best_alpha_C) * gaussian_sphere_density(sigma_halo)
im6 = ax6.pcolormesh(rho_arr, z_arr, rho_mix.T, cmap='inferno', shading='auto')
ax6.set_xlabel('rho [fm]')
ax6.set_ylabel('z [fm]')
ax6.set_title(f'Torus + halo (alpha={best_alpha_C:.2f})')
ax6.set_aspect('equal')
plt.colorbar(im6, ax=ax6, shrink=0.7, label='charge density')

fig.suptitle('Can Quantum Averaging Hide a Toroidal Proton?', fontsize=15, y=1.01)
plt.tight_layout()
plt.savefig('toroidal_proton_quantum_averaged.png', dpi=150, bbox_inches='tight')
print(f"Plot saved to toroidal_proton_quantum_averaged.png")
plt.show()
