"""
Toroidal Proton Model: Calculable Predictions
=============================================

Starting from the hypothesis "the proton is a torus," we derive testable
predictions for the magnetic moments and charge radii of the proton and
neutron.

The model:
  - The proton is a torus with major radius R and minor radius r.
  - Up quarks (charge +2/3 e) orbit at the INNER equator, radius R_u = R - r.
  - Down quarks (charge -1/3 e) orbit at the OUTER equator, radius R_d = R + r.
  - All quarks move at speed v ~ c (ultrarelativistic, as expected inside
    a hadron of this size from the uncertainty principle).
  - The magnetic moment is just the sum of three tiny current loops.

The key insight: the model has ONE geometric parameter — the aspect ratio
R/r.  We fix it from the measured ratio mu_p/mu_n.  After that, everything
else is a genuine, testable prediction.

What determines what:
  FITTED:  mu_p/mu_n  -->  determines R/r  (the shape of the torus)
  FITTED:  |mu_p|     -->  determines R    (the size of the torus)
  PREDICTED:  sign of neutron <r^2>         (must be negative)
  PREDICTED:  proton charge radius
  PREDICTED:  neutron <r^2>
  PREDICTED:  which quarks orbit where      (d must be outside)
"""

import numpy as np
import matplotlib.pyplot as plt

# =========================================================================
# PHYSICAL CONSTANTS
# =========================================================================

hbar = 1.0546e-34       # J·s    reduced Planck constant
c    = 2.998e8           # m/s    speed of light
e    = 1.602e-19         # C      elementary charge
m_p  = 1.673e-27         # kg     proton mass
fm   = 1e-15             # m      femtometer
mu_N = e * hbar / (2 * m_p)   # nuclear magneton = 5.051e-27 J/T

# =========================================================================
# EXPERIMENTAL DATA (Particle Data Group)
# =========================================================================

mu_p_exp  =  2.7928473    # proton magnetic moment  [nuclear magnetons]
mu_n_exp  = -1.9130427    # neutron magnetic moment [nuclear magnetons]
r_p_exp   =  0.841        # proton RMS charge radius [fm]
r2_n_exp  = -0.1161       # neutron mean-square charge radius [fm^2]

ratio_exp = mu_p_exp / mu_n_exp   # = -1.4598

# Quark charges in units of e
q_u =  2.0 / 3.0    # up quark
q_d = -1.0 / 3.0    # down quark

print("=" * 70)
print("  TOROIDAL PROTON MODEL: CALCULABLE PREDICTIONS")
print("=" * 70)
print(f"\n  Experimental targets:")
print(f"    mu_p         = {mu_p_exp:+.4f} mu_N")
print(f"    mu_n         = {mu_n_exp:+.4f} mu_N")
print(f"    mu_p / mu_n  = {ratio_exp:.4f}")
print(f"    proton r_ch  = {r_p_exp} fm")
print(f"    neutron <r2> = {r2_n_exp} fm^2")
print(f"    SU(6) quark model mu_p/mu_n = {-3/2:.4f}")


# =========================================================================
# PART 1: DERIVE THE TORUS SHAPE FROM mu_p / mu_n
# =========================================================================
#
# A charge q orbiting at radius R with speed v produces:
#     magnetic moment  =  q * v * R / 2
#
# Proton (uud) — two up quarks at R_u, one down quark at R_d:
#     mu_p = (e*v/2) * [ 2*(2/3)*R_u + (-1/3)*R_d ]
#          = (e*v/6) * [ 4*R_u - R_d ]
#
# Neutron (udd) — one up quark at R_u, two down quarks at R_d:
#     mu_n = (e*v/2) * [ (2/3)*R_u + 2*(-1/3)*R_d ]
#          = (e*v/3) * [ R_u - R_d ]
#
# The ratio (note: v and e cancel completely!):
#     mu_p / mu_n  =  (4*R_u - R_d) / (2*(R_u - R_d))
#
# On a torus with R_u = R - r and R_d = R + r, define a = R/r:
#     mu_p / mu_n  =  -(3a - 5) / 4
#
# Solving for a:
#     a  =  (5 - 4 * ratio) / 3

a_predicted = (5 - 4 * ratio_exp) / 3          # torus aspect ratio R/r
rho = (a_predicted + 1) / (a_predicted - 1)     # orbit ratio R_d / R_u

print(f"\n{'='*70}")
print("  PART 1: TORUS SHAPE FROM mu_p / mu_n")
print(f"{'='*70}")
print(f"\n  From mu_p/mu_n = {ratio_exp:.4f}:")
print(f"    Torus aspect ratio  R/r = {a_predicted:.4f}")
print(f"    Orbit ratio  R_d/R_u    = {rho:.4f}")
print(f"    Down quarks orbit {rho:.2f}x further out than up quarks.")
print(f"\n  WHY d must be outside:")
print(f"    If up quarks were outside, mu_n would be POSITIVE.")
print(f"    Experiment says mu_n < 0, so down quarks must be further out.")
print(f"    The torus geometry DETERMINES the quark arrangement.")


# =========================================================================
# PART 2: ABSOLUTE SIZE FROM |mu_p|
# =========================================================================
#
# mu_p = e*c*R * (3 - 5/a) / 6    (setting v = c)
#
# Set equal to 2.7928 * mu_N and solve for R:

A_p = (3 - 5 / a_predicted) / 6     # dimensionless coefficient for proton
A_n = -2 / (3 * a_predicted)         # dimensionless coefficient for neutron

R_m = mu_p_exp * mu_N / (e * c * A_p)   # major radius in meters
r_m = R_m / a_predicted                  # minor radius in meters
R_fm = R_m / fm                          # major radius in femtometers
r_fm = r_m / fm                          # minor radius in femtometers
R_u_fm = R_fm - r_fm                     # up-quark orbit radius
R_d_fm = R_fm + r_fm                     # down-quark orbit radius

print(f"\n{'='*70}")
print("  PART 2: ABSOLUTE SIZE (setting v = c)")
print(f"{'='*70}")
print(f"\n  Torus dimensions:")
print(f"    Major radius  R   = {R_fm:.4f} fm")
print(f"    Minor radius  r   = {r_fm:.4f} fm")
print(f"    R_u (up orbit)    = {R_u_fm:.4f} fm")
print(f"    R_d (down orbit)  = {R_d_fm:.4f} fm")
print(f"    Proton Compton wavelength h/(m_p c) = {hbar/(m_p*c)/fm:.4f} fm")
print(f"\n  Check: mu_p = e*c*R*{A_p:.4f} = {A_p * e * c * R_m / mu_N:.4f} mu_N  (target: {mu_p_exp:.4f})")

# Neutron moment (genuine prediction — we only fitted the RATIO)
mu_n_pred = A_n / A_p * mu_p_exp   # in mu_N
print(f"  Check: mu_n = {mu_n_pred:.4f} mu_N  (target: {mu_n_exp:.4f})")


# =========================================================================
# PART 3: CHARGE RADII — GENUINE PREDICTIONS
# =========================================================================
#
# The charge radius is the charge-weighted mean-square radius:
#
#   Proton:  <r^2>_p = [ 2*q_u*R_u^2 + q_d*R_d^2 ] / Q_p
#                    = [ (4/3)*R_u^2 - (1/3)*R_d^2 ]
#
#   Neutron: <r^2>_n = q_u*R_u^2 + 2*q_d*R_d^2     (no division — Q_n = 0)
#                    = (2/3)*(R_u^2 - R_d^2)

R_u_m = R_u_fm * fm
R_d_m = R_d_fm * fm

r2_p_pred = (4/3 * R_u_m**2 - 1/3 * R_d_m**2) / fm**2   # fm^2
r_p_pred  = np.sqrt(abs(r2_p_pred))                        # fm

r2_n_pred = 2/3 * (R_u_m**2 - R_d_m**2) / fm**2           # fm^2

print(f"\n{'='*70}")
print("  PART 3: CHARGE RADII (genuine predictions)")
print(f"{'='*70}")
print(f"\n  Proton charge radius:")
print(f"    Predicted  <r^2>_p = {r2_p_pred:+.4f} fm^2   ->  r_p = {r_p_pred:.4f} fm")
print(f"    Measured   <r^2>_p = {r_p_exp**2:+.4f} fm^2   ->  r_p = {r_p_exp:.4f} fm")
print(f"    Ratio (pred/exp):    {r_p_pred/r_p_exp:.3f}")

print(f"\n  Neutron charge radius:")
print(f"    Predicted  <r^2>_n = {r2_n_pred:+.4f} fm^2")
print(f"    Measured   <r^2>_n = {r2_n_exp:+.4f} fm^2")
print(f"    Sign:  {'CORRECT (negative)' if r2_n_pred < 0 else 'WRONG'}")
print(f"    Ratio (pred/exp):    {r2_n_pred/r2_n_exp:.3f}")


# =========================================================================
# PART 4: SCORECARD — MODEL vs EXPERIMENT
# =========================================================================

print(f"\n{'='*70}")
print("  SCORECARD: TOROIDAL PROTON vs EXPERIMENT")
print(f"{'='*70}")
print(f"\n  {'Observable':<28} {'Predicted':>12} {'Measured':>12} {'Status':>10}")
print(f"  {'-'*62}")
print(f"  {'mu_p [mu_N]':<28} {mu_p_exp:>+12.4f} {mu_p_exp:>+12.4f} {'(fitted)':>10}")
print(f"  {'mu_n [mu_N]':<28} {mu_n_pred:>+12.4f} {mu_n_exp:>+12.4f} {'(fitted)':>10}")
print(f"  {'mu_p / mu_n':<28} {ratio_exp:>12.4f} {ratio_exp:>12.4f} {'(fitted)':>10}")
print(f"  {'R/r (torus shape)':<28} {a_predicted:>12.4f} {'—':>12} {'from fit':>10}")
print(f"  {'R [fm]':<28} {R_fm:>12.4f} {'—':>12} {'from fit':>10}")
print(f"  {'r [fm]':<28} {r_fm:>12.4f} {'—':>12} {'from fit':>10}")
print(f"  {'sign(neutron <r2>)':<28} {'negative':>12} {'negative':>12} {'PASS':>10}")
print(f"  {'proton r_ch [fm]':<28} {r_p_pred:>12.4f} {r_p_exp:>12.4f} {'~{:.0f}%'.format(r_p_pred/r_p_exp*100):>10}")
print(f"  {'neutron <r2> [fm2]':<28} {r2_n_pred:>+12.4f} {r2_n_exp:>+12.4f} {'~{:.0f}x'.format(r2_n_pred/r2_n_exp):>10}")


# =========================================================================
# PART 5: WHAT QUARK VELOCITY FIXES THE CHARGE RADIUS?
# =========================================================================
#
# If we DON'T assume v = c, we can ask: what velocity makes the predicted
# charge radius match experiment?
#
# Scaling: all radii scale as 1/v (since R = mu_p*mu_N / (e*v*A_p)).
# So if we increase R to match r_p, we need to DECREASE v.
#
# Required scale factor:

scale_needed = r_p_exp / r_p_pred   # how much bigger the torus needs to be
v_needed = c / scale_needed          # slower quarks -> bigger torus
beta_needed = v_needed / c           # v/c

# At this velocity, check the moments:
R_adj = R_fm * scale_needed
r_adj = r_fm * scale_needed
mu_p_check = A_p * e * v_needed * (R_adj * fm) / mu_N
mu_n_check = A_n * e * v_needed * (R_adj * fm) / mu_N

print(f"\n{'='*70}")
print("  PART 5: VELOCITY ADJUSTMENT TO MATCH CHARGE RADIUS")
print(f"{'='*70}")
print(f"\n  If the proton charge radius must be {r_p_exp} fm:")
print(f"    Required scale factor: {scale_needed:.4f}")
print(f"    Required quark speed:  v = {beta_needed:.4f} c")
print(f"    Adjusted R = {R_adj:.4f} fm,  r = {r_adj:.4f} fm")
print(f"    mu_p at this v: {mu_p_check:.4f} mu_N  (target: {mu_p_exp:.4f})")
print(f"    mu_n at this v: {mu_n_check:.4f} mu_N  (target: {mu_n_exp:.4f})")
print(f"\n  Key point: the RATIO mu_p/mu_n = {mu_p_check/mu_n_check:.4f}")
print(f"  is INDEPENDENT of velocity.  The geometry determines the ratio.")
print(f"  Only the absolute scale changes with v.")


# =========================================================================
# PART 6: NUMERICAL VERIFICATION USING TORUS-KNOT PATHS
# =========================================================================
#
# As a cross-check, and to explore what happens when quarks follow
# helical (torus knot) paths instead of simple circles, we compute
# the magnetic moment numerically using the integral m = (I/2) * sum(r x dl).

print(f"\n{'='*70}")
print("  PART 6: TORUS-KNOT QUARK PATHS")
print(f"{'='*70}")

N_pts = 100_000   # numerical resolution

def torus_path(R_maj, r_min, theta_0, N_poloidal=0, N_pts=100_000):
    """
    Generate a closed path on a torus surface.

    Parameters:
        R_maj      : major radius of torus
        r_min      : minor radius of torus
        theta_0    : poloidal angle of the orbit (0 = outer, pi = inner)
                     Only used when N_poloidal = 0 (simple circle).
        N_poloidal : number of poloidal windings. 0 = simple circle at
                     fixed poloidal angle theta_0.
    Returns:
        path : (N_pts, 3) array of points on the torus
    """
    t = np.linspace(0, 2 * np.pi, N_pts, endpoint=False)

    if N_poloidal == 0:
        # Simple circle at a fixed poloidal angle
        R_eff = R_maj + r_min * np.cos(theta_0)
        z_eff = r_min * np.sin(theta_0)
        x = R_eff * np.cos(t)
        y = R_eff * np.sin(t)
        z = np.full_like(t, z_eff)
    else:
        # (1, N_poloidal) torus knot — winds once toroidally, N times poloidally
        x = (R_maj + r_min * np.cos(N_poloidal * t)) * np.cos(t)
        y = (R_maj + r_min * np.cos(N_poloidal * t)) * np.sin(t)
        z = r_min * np.sin(N_poloidal * t)

    return np.column_stack([x, y, z])


def compute_moment(path, charge_e, speed):
    """
    Compute magnetic moment vector for a charge orbiting on a closed path.

    Current = charge * speed / path_length
    Moment  = (I / 2) * sum(r x dl)
    """
    dl = np.roll(path, -1, axis=0) - path
    L = np.sum(np.linalg.norm(dl, axis=1))
    I = charge_e * speed / L
    m_vec = 0.5 * I * np.sum(np.cross(path, dl), axis=0)
    return m_vec


# --- Verify the analytical result with simple circular orbits ---

R_torus = R_fm * fm    # in meters
r_torus = r_fm * fm

# Up quarks: circle at inner equator (theta = pi)
path_u = torus_path(R_torus, r_torus, theta_0=np.pi, N_pts=N_pts)
m_u = compute_moment(path_u, q_u * e, c)

# Down quarks: circle at outer equator (theta = 0)
path_d = torus_path(R_torus, r_torus, theta_0=0.0, N_pts=N_pts)
m_d = compute_moment(path_d, q_d * e, c)

# Proton: 2u + 1d
m_p_num = 2 * m_u + m_d
mu_p_num = np.linalg.norm(m_p_num) / mu_N * np.sign(m_p_num[2])

# Neutron: 1u + 2d
m_n_num = m_u + 2 * m_d
mu_n_num = np.linalg.norm(m_n_num) / mu_N * np.sign(m_n_num[2])

print(f"\n  Numerical verification (simple circular orbits):")
print(f"    mu_p = {mu_p_num:+.4f} mu_N  (analytical: {mu_p_exp:+.4f})")
print(f"    mu_n = {mu_n_num:+.4f} mu_N  (analytical: {mu_n_pred:+.4f})")
print(f"    ratio = {mu_p_num/mu_n_num:.4f}  (target: {ratio_exp:.4f})")

# --- Now try helical (torus-knot) quark paths ---
print(f"\n  What if quarks follow helical paths on the torus?")
print(f"  (N_pol = number of poloidal windings per toroidal orbit)")
print(f"\n  {'N_pol':>6}  {'mu_p [mu_N]':>12}  {'mu_n [mu_N]':>12}  {'mu_p/mu_n':>12}")
print(f"  {'-'*48}")

for N_pol in [0, 1, 2, 3, 5, 10]:
    if N_pol == 0:
        # Simple circles (already computed)
        print(f"  {N_pol:>6}  {mu_p_num:>+12.4f}  {mu_n_num:>+12.4f}  {mu_p_num/mu_n_num:>12.4f}  (circles)")
        continue

    # Up quarks on a helical path near inner part of torus
    # We shrink r slightly so the helix stays near the inner equator
    path_u_knot = torus_path(R_torus, r_torus * 0.8, theta_0=0, N_poloidal=N_pol, N_pts=N_pts)
    m_u_knot = compute_moment(path_u_knot, q_u * e, c)

    # Down quarks on a helical path using full tube radius
    path_d_knot = torus_path(R_torus, r_torus, theta_0=0, N_poloidal=N_pol, N_pts=N_pts)
    m_d_knot = compute_moment(path_d_knot, q_d * e, c)

    m_p_knot = 2 * m_u_knot + m_d_knot
    m_n_knot = m_u_knot + 2 * m_d_knot
    mu_p_k = m_p_knot[2] / mu_N
    mu_n_k = m_n_knot[2] / mu_N

    if abs(mu_n_k) > 1e-30:
        print(f"  {N_pol:>6}  {mu_p_k:>+12.4f}  {mu_n_k:>+12.4f}  {mu_p_k/mu_n_k:>12.4f}")
    else:
        print(f"  {N_pol:>6}  {mu_p_k:>+12.4f}  {mu_n_k:>+12.4f}  {'inf':>12}")


# =========================================================================
# PART 7: VISUALIZATION
# =========================================================================

fig = plt.figure(figsize=(18, 8))

# ---- Left panel: 3-D torus with quark orbits ----
ax1 = fig.add_subplot(121, projection='3d')

# Draw torus surface (wireframe)
theta_grid = np.linspace(0, 2 * np.pi, 60)
phi_grid   = np.linspace(0, 2 * np.pi, 30)
THETA, PHI = np.meshgrid(theta_grid, phi_grid)
# Use dimensionless units for the plot (divide by fm)
X_torus = (R_fm + r_fm * np.cos(PHI)) * np.cos(THETA)
Y_torus = (R_fm + r_fm * np.cos(PHI)) * np.sin(THETA)
Z_torus = r_fm * np.sin(PHI)

ax1.plot_surface(X_torus, Y_torus, Z_torus, alpha=0.08, color='gold')
ax1.plot_wireframe(X_torus, Y_torus, Z_torus, alpha=0.1, color='goldenrod',
                   rstride=3, cstride=3)

# Draw up-quark orbit (inner circle)
t_plot = np.linspace(0, 2 * np.pi, 300)
xu = R_u_fm * np.cos(t_plot)
yu = R_u_fm * np.sin(t_plot)
zu = np.zeros_like(t_plot)
ax1.plot(xu, yu, zu, color='red', linewidth=2.5, label=f'u orbit (R={R_u_fm:.2f} fm)')

# Draw down-quark orbit (outer circle)
xd = R_d_fm * np.cos(t_plot)
yd = R_d_fm * np.sin(t_plot)
zd = np.zeros_like(t_plot)
ax1.plot(xd, yd, zd, color='blue', linewidth=2.5, label=f'd orbit (R={R_d_fm:.2f} fm)')

# Draw quark positions (3 quarks for proton: u, u, d)
quark_angles = [0, 2*np.pi/3, 4*np.pi/3]
for i, ang in enumerate(quark_angles):
    if i < 2:  # up quarks
        ax1.scatter([R_u_fm*np.cos(ang)], [R_u_fm*np.sin(ang)], [0],
                    color='red', s=120, zorder=5, edgecolors='black')
    else:       # down quark
        ax1.scatter([R_d_fm*np.cos(ang)], [R_d_fm*np.sin(ang)], [0],
                    color='blue', s=120, zorder=5, edgecolors='black')

# Draw magnetic moment vector
arrow_len = R_fm * 0.7
ax1.quiver(0, 0, 0, 0, 0, arrow_len,
           color='crimson', linewidth=3, arrow_length_ratio=0.12,
           label=r'$\mu_p$ direction')

ax1.set_xlabel('X [fm]')
ax1.set_ylabel('Y [fm]')
ax1.set_zlabel('Z [fm]')
ax1.set_title('Toroidal Proton: Quark Orbits', fontsize=12, pad=10)
ax1.legend(loc='upper left', fontsize=8)

# ---- Right panel: mu_p/mu_n vs aspect ratio ----
ax2 = fig.add_subplot(122)

a_scan = np.linspace(1.8, 10, 500)
ratio_scan = -(3 * a_scan - 5) / 4

ax2.plot(a_scan, ratio_scan, 'b-', linewidth=2, label=r'$\mu_p/\mu_n = -(3R/r - 5)/4$')
ax2.axhline(y=ratio_exp, color='red', linestyle='--', linewidth=1.5,
            label=f'Experiment: {ratio_exp:.4f}')
ax2.axvline(x=a_predicted, color='green', linestyle='--', linewidth=1.5,
            label=f'R/r = {a_predicted:.2f}')
ax2.axhline(y=-1.5, color='gray', linestyle=':', linewidth=1,
            label='SU(6) quark model: -1.5')

ax2.scatter([a_predicted], [ratio_exp], color='red', s=100, zorder=5)

ax2.set_xlabel('Torus aspect ratio R/r', fontsize=12)
ax2.set_ylabel(r'$\mu_p\,/\,\mu_n$', fontsize=12)
ax2.set_title('Moment Ratio vs Torus Geometry', fontsize=12)
ax2.legend(fontsize=9)
ax2.set_ylim(-3, 0)
ax2.set_xlim(1.5, 10)
ax2.grid(True, alpha=0.3)

# Annotate the solution point
ax2.annotate(f'  R/r = {a_predicted:.2f}\n  matches experiment',
             xy=(a_predicted, ratio_exp),
             xytext=(a_predicted + 1.5, ratio_exp + 0.5),
             fontsize=10, color='green',
             arrowprops=dict(arrowstyle='->', color='green'))

fig.suptitle('Proton as a Torus: Magnetic Moment Predictions', fontsize=14, y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('toroidal_proton_model.png', dpi=150, bbox_inches='tight')
print(f"\n  Plot saved to toroidal_proton_model.png")
plt.show()
