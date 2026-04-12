"""
Magnetic Moment of a Trefoil Knot Current Loop
================================================

This script calculates the magnetic dipole moment of a current-carrying wire
shaped into a trefoil knot, then compares it to:
  (a) a simple circular loop of the same total wire length, and
  (b) a toroidal helical coil (wire wound around a donut shape) of the same
      total wire length.

All three carry the same current I = 1.

Background
----------
The magnetic dipole moment of a current loop is:

    m = (I / 2) * integral( r x dl )

where r is the position vector along the wire and dl is an infinitesimal
tangent element.  For a planar loop this simplifies to m = I * A (current
times enclosed area), but a trefoil knot is a genuinely 3-D curve, so we
must evaluate the vector integral numerically.

A trefoil knot can be parameterised as:
    x(t) = sin(t) + 2 sin(2t)
    y(t) = cos(t) - 2 cos(2t)
    z(t) = -sin(3t)
with t in [0, 2*pi).

A toroidal helix (wire wound around a torus) is parameterised as:
    x(t) = (R + r cos(N_p t)) cos(t)
    y(t) = (R + r cos(N_p t)) sin(t)
    z(t) = r sin(N_p t)
where R is the major radius, r is the minor radius, and N_p is the number
of poloidal (small-circle) windings the wire makes as it goes once around
the torus.  The interesting physics: as N_p grows the poloidal loops' moments
increasingly cancel, so the net moment drops toward zero -- the hallmark of
an ideal toroidal solenoid whose field is confined inside the donut.
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Step 1: Define the trefoil knot as a parametric 3-D curve
# ---------------------------------------------------------------------------
# We sample the curve at a large number of points so the numerical
# integration is accurate.  More points = more accuracy but slower.

N = 100_000                        # number of sample points along the knot
t = np.linspace(0, 2 * np.pi, N, endpoint=False)  # parameter values

# Standard trefoil knot parametrisation
x = np.sin(t) + 2 * np.sin(2 * t)
y = np.cos(t) - 2 * np.cos(2 * t)
z = -np.sin(3 * t)

# Stack into an (N, 3) array: each row is one point [x, y, z] on the knot.
r = np.column_stack([x, y, z])     # position vectors along the curve

# ---------------------------------------------------------------------------
# Step 2: Compute the tangent (dl) vectors along the curve
# ---------------------------------------------------------------------------
# dl[i] is the tiny vector from point i to point i+1.  We use np.roll to
# wrap around so the curve closes back on itself (it's a closed knot).

dl = np.roll(r, -1, axis=0) - r    # dl[i] = r[i+1] - r[i]

# ---------------------------------------------------------------------------
# Step 3: Compute the total arc length of the trefoil knot
# ---------------------------------------------------------------------------
# The arc length is the sum of all the tiny step lengths |dl|.

segment_lengths = np.linalg.norm(dl, axis=1)   # length of each tiny step
L_trefoil = np.sum(segment_lengths)             # total wire length

print(f"Total arc length of the trefoil knot: {L_trefoil:.4f} (dimensionless units)")

# ---------------------------------------------------------------------------
# Step 4: Calculate the magnetic dipole moment of the trefoil knot
# ---------------------------------------------------------------------------
# The magnetic moment of a current loop carrying current I is:
#
#     m_vec = (I / 2) * SUM_over_curve( r_i x dl_i )
#
# where "x" is the cross product.  We set I = 1 (unit current) so the
# result is purely geometric -- it depends only on the shape of the knot.

I_current = 1.0                     # current (arbitrary units)

# Cross product of each position vector with its tangent element
cross_products = np.cross(r, dl)    # shape (N, 3)

# Sum all the tiny contributions and multiply by I/2
m_trefoil_vec = 0.5 * I_current * np.sum(cross_products, axis=0)

# The magnitude (scalar) of the magnetic moment
m_trefoil_mag = np.linalg.norm(m_trefoil_vec)

print(f"\nTrefoil knot magnetic moment vector: [{m_trefoil_vec[0]:+.6f}, "
      f"{m_trefoil_vec[1]:+.6f}, {m_trefoil_vec[2]:+.6f}]")
print(f"Trefoil knot magnetic moment magnitude: {m_trefoil_mag:.6f}")

# ---------------------------------------------------------------------------
# Step 5: Calculate the magnetic moment of a circular loop with the SAME
#          total wire length and the SAME current
# ---------------------------------------------------------------------------
# A circle of circumference L has radius R = L / (2*pi).
# Its magnetic moment is simply  m = I * pi * R^2  (current times area).

R_circle = L_trefoil / (2 * np.pi)              # radius of the equivalent circle
A_circle = np.pi * R_circle**2                   # area of the circle
m_circle = I_current * A_circle                  # magnetic moment magnitude

print(f"\nEquivalent circular loop radius: {R_circle:.4f}")
print(f"Circular loop magnetic moment magnitude: {m_circle:.6f}")

# ---------------------------------------------------------------------------
# Step 6: Calculate the magnetic moment of a TOROIDAL HELIX (wire wound
#          around a donut) with the SAME total wire length and current
# ---------------------------------------------------------------------------
# A toroidal helix is a wire that spirals around a torus, making N_p
# poloidal (small-circle) turns as it completes one full toroidal (big-circle)
# turn.  We choose an aspect ratio R/r = 3 (a "fat" donut), then uniformly
# scale the whole torus so the total wire length matches the trefoil.
#
# Physics insight: each poloidal loop acts like a small current loop whose
# magnetic moment points in the local toroidal (tangential) direction.
# Because these small moments point in different directions all around
# the ring, they partially cancel.  With more poloidal turns the cancellation
# becomes more complete -- an ideal toroidal solenoid (infinite N_p) has
# ZERO net dipole moment because its field is entirely trapped inside.

def compute_torus_moment(N_poloidal, L_target, N_pts=100_000):
    """Compute the magnetic moment of a toroidal helix with N_poloidal
    poloidal windings, scaled so total arc length equals L_target."""

    # Start with an unscaled torus: aspect ratio R0/r0 = 3
    R0 = 3.0    # major radius (centre of tube to centre of torus)
    r0 = 1.0    # minor radius (radius of the tube cross-section)

    t = np.linspace(0, 2 * np.pi, N_pts, endpoint=False)

    # Parametric curve on the torus surface
    x_t = (R0 + r0 * np.cos(N_poloidal * t)) * np.cos(t)
    y_t = (R0 + r0 * np.cos(N_poloidal * t)) * np.sin(t)
    z_t = r0 * np.sin(N_poloidal * t)

    r_t = np.column_stack([x_t, y_t, z_t])
    dl_t = np.roll(r_t, -1, axis=0) - r_t
    L_raw = np.sum(np.linalg.norm(dl_t, axis=1))

    # Uniformly scale so total path length matches L_target
    scale = L_target / L_raw
    r_t *= scale
    dl_t *= scale

    # Magnetic moment: m = (I/2) * sum(r x dl)
    m_vec = 0.5 * I_current * np.sum(np.cross(r_t, dl_t), axis=0)
    return m_vec, r_t

# We'll compute the torus moment for several winding numbers to show
# how the moment changes as the coil wraps more tightly.
poloidal_counts = [1, 2, 5, 10, 20, 50]

print("\n--- Toroidal helix comparison (same wire length, same current) ---")
print(f"{'N_poloidal':>10}  {'|m_torus|':>12}  {'m_torus/m_circle':>16}  {'m_torus/m_trefoil':>18}")
print("-" * 65)

torus_moments = {}
for Np in poloidal_counts:
    m_torus_vec, _ = compute_torus_moment(Np, L_trefoil)
    m_torus_mag = np.linalg.norm(m_torus_vec)
    torus_moments[Np] = (m_torus_vec, m_torus_mag)
    print(f"{Np:>10}  {m_torus_mag:>12.6f}  {m_torus_mag/m_circle:>16.6f}  {m_torus_mag/m_trefoil_mag:>18.6f}")

# For the 3-D plot we'll use N_p = 10 as a representative toroidal coil.
N_p_plot = 10
m_torus_plot_vec, r_torus_plot = compute_torus_moment(N_p_plot, L_trefoil)
m_torus_plot_mag = np.linalg.norm(m_torus_plot_vec)

# ---------------------------------------------------------------------------
# Step 7: Compare all three moments
# ---------------------------------------------------------------------------
ratio = m_trefoil_mag / m_circle
ratio_torus = m_torus_plot_mag / m_circle

print(f"\n--- Summary ---")
print(f"Circular loop moment:           {m_circle:.6f}")
print(f"Trefoil knot moment:            {m_trefoil_mag:.6f}  ({ratio*100:.1f}% of circle)")
print(f"Toroidal helix (Np={N_p_plot}) moment: {m_torus_plot_mag:.6f}  ({ratio_torus*100:.1f}% of circle)")
print(f"\nRatio  m_trefoil / m_circle = {ratio:.6f}")
print(f"Ratio  m_torus   / m_circle = {ratio_torus:.6f}")
print(f"Ratio  m_trefoil / m_torus  = {m_trefoil_mag/m_torus_plot_mag:.6f}")
print(f"\nKey insight: the toroidal solenoid's poloidal loops cancel each other's")
print(f"moments.  As N_poloidal -> infinity, the net moment -> 0.  The trefoil,")
print(f"being a (2,3) torus knot, sits between a simple circle and a tightly")
print(f"wound toroid.")

# ---------------------------------------------------------------------------
# Step 8: Plot both curves in 3-D with their magnetic moment vectors
# ---------------------------------------------------------------------------

fig = plt.figure(figsize=(18, 8))

# --- Left panel: Trefoil knot ---
ax1 = fig.add_subplot(121, projection='3d')

ax1.plot(x, y, z, color='royalblue', linewidth=1.2, label='Trefoil knot')

# Magnetic moment arrow (scaled for visibility)
arrow_scale = 3.0 / m_trefoil_mag
mx, my, mz = m_trefoil_vec * arrow_scale
ax1.quiver(0, 0, 0, mx, my, mz,
           color='crimson', linewidth=3, arrow_length_ratio=0.12,
           label='Magnetic moment')
ax1.scatter([0], [0], [0], color='black', s=30, zorder=5)

ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.set_title('Trefoil Knot (2,3 torus knot)', fontsize=12, pad=10)
ax1.legend(loc='upper left', fontsize=9)

# --- Right panel: Toroidal helix ---
ax2 = fig.add_subplot(122, projection='3d')

xt, yt, zt = r_torus_plot[:, 0], r_torus_plot[:, 1], r_torus_plot[:, 2]
ax2.plot(xt, yt, zt, color='seagreen', linewidth=1.0,
         label=f'Toroidal helix (Np={N_p_plot})')

# Magnetic moment arrow for the torus
if m_torus_plot_mag > 1e-10:
    arrow_scale_t = 3.0 / m_torus_plot_mag
else:
    arrow_scale_t = 1.0
mtx, mty, mtz = m_torus_plot_vec * arrow_scale_t
ax2.quiver(0, 0, 0, mtx, mty, mtz,
           color='crimson', linewidth=3, arrow_length_ratio=0.12,
           label='Magnetic moment')
ax2.scatter([0], [0], [0], color='black', s=30, zorder=5)

ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.set_title(f'Toroidal Helix ({N_p_plot} poloidal turns)', fontsize=12, pad=10)
ax2.legend(loc='upper left', fontsize=9)

# --- Shared annotation ---
textstr = (
    f"Same wire length L = {L_trefoil:.2f},  same current I = {I_current}\n"
    f"|m_circle|  = {m_circle:.2f}   (reference: best possible)\n"
    f"|m_trefoil| = {m_trefoil_mag:.2f}   ({ratio*100:.1f}% of circle)\n"
    f"|m_torus|   = {m_torus_plot_mag:.2f}   ({ratio_torus*100:.1f}% of circle)"
)
fig.text(0.5, 0.01, textstr, fontsize=11, fontfamily='monospace',
         ha='center', va='bottom',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

fig.suptitle('Magnetic Dipole Moment: Trefoil Knot vs Toroidal Helix',
             fontsize=15, y=0.98)
plt.tight_layout(rect=[0, 0.08, 1, 0.95])
plt.savefig('trefoil_magnetic_moment.png', dpi=150, bbox_inches='tight')
print("\nPlot saved to trefoil_magnetic_moment.png")
plt.show()
