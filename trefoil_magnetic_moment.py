"""
Magnetic Moment of a Trefoil Knot Current Loop
================================================

This script calculates the magnetic dipole moment of a current-carrying wire
shaped into a trefoil knot, then compares it to a simple circular loop of the
same total wire length carrying the same current.

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
# Step 6: Compare the two moments
# ---------------------------------------------------------------------------
ratio = m_trefoil_mag / m_circle

print(f"\nRatio  m_trefoil / m_circle = {ratio:.6f}")
print(f"The trefoil knot produces only ~{ratio*100:.1f}% of the magnetic "
      f"moment of a circular loop of the same wire length.")

# ---------------------------------------------------------------------------
# Step 7: Plot the trefoil knot in 3-D with the magnetic moment vector
# ---------------------------------------------------------------------------

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Draw the knot itself, coloured by the parameter t so you can see how
# the curve winds around.
ax.plot(x, y, z, color='royalblue', linewidth=1.2, label='Trefoil knot')

# Draw the magnetic moment vector as a thick arrow starting at the origin.
# We scale the arrow so it is visually comparable to the knot size.
arrow_scale = 3.0 / m_trefoil_mag   # normalise then scale to length 3
mx, my, mz = m_trefoil_vec * arrow_scale
ax.quiver(0, 0, 0, mx, my, mz,
          color='crimson', linewidth=3, arrow_length_ratio=0.12,
          label=f'Magnetic moment direction')

# Mark the origin
ax.scatter([0], [0], [0], color='black', s=30, zorder=5)

# Labels and legend
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Trefoil Knot Current Loop & Magnetic Dipole Moment',
             fontsize=13, pad=15)
ax.legend(loc='upper left', fontsize=10)

# Add a text box with the key result
textstr = (f"|m_trefoil| / |m_circle| = {ratio:.4f}\n"
           f"m direction = ({m_trefoil_vec[0]/m_trefoil_mag:+.3f}, "
           f"{m_trefoil_vec[1]/m_trefoil_mag:+.3f}, "
           f"{m_trefoil_vec[2]/m_trefoil_mag:+.3f})")
fig.text(0.02, 0.02, textstr, fontsize=10, fontfamily='monospace',
         verticalalignment='bottom',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('trefoil_magnetic_moment.png', dpi=150, bbox_inches='tight')
print("\nPlot saved to trefoil_magnetic_moment.png")
plt.show()
