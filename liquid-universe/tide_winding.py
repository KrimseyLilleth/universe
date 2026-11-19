"""
Count how the tide's timing swirls around every amphidromic point AND every landmass.

What it checks
--------------
At an amphidromic point, the time of high tide turns through a full cycle as you walk
around it. A landmass can do the same thing: walk once around its coast and high tide's
timing can lap a whole number of times, with no still point in the water. That is a
"hole" doing the job of an amphidromic point, the same way the hole of a torus can.

Topology says the swirls in the water plus the swirls around the coasts must add up to
exactly zero over the whole globe. This script measures both and prints the totals.

Run:  python3 tide_winding.py            (needs numpy and scipy)
Data: amphidromic_got410c_data.js (NASA GOT4.10c, 0.5 degree grid; rows lat -90..90,
      columns lon -180..179.5). +1 = counter-clockwise seen from above, -1 = clockwise.
"""
import base64
import json
import os

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
PLACES = {  # a point inside each landmass we want to name
    "Africa + Eurasia": (5, 20), "The Americas": (40, -100), "Antarctica": (-85, 0),
    "Greenland": (75, -40), "Australia": (-25, 135), "New Zealand (South Is.)": (-44, 170.5),
    "Madagascar": (-20, 46.5), "Japan (Honshu)": (36, 138), "Iceland": (65, -18),
    "Ireland": (53, -8), "Borneo": (0, 114), "Sri Lanka": (7.5, 80.7),
}


def load():
    s = open(os.path.join(HERE, "amphidromic_got410c_data.js")).read()
    return json.loads(s[s.index("{"): s.rindex("}") + 1])


def wrap(x):
    return np.angle(np.exp(1j * x))  # into -pi..pi


def main():
    data = load()
    for name, c in data["constituents"].items():
        f = c["field"]
        nlat, nlon = f["nlat"], f["nlon"]
        raw = np.frombuffer(base64.b64decode(f["b64"]), dtype="<i2").reshape(nlat, nlon, 2)
        land = raw[..., 0] == -32768
        phase = np.arctan2(raw[..., 1].astype(float), raw[..., 0].astype(float))

        # phase step along each grid edge (set to 0 if the edge touches land)
        east = wrap(np.roll(phase, -1, 1) - phase)
        east[land | np.roll(land, -1, 1)] = 0
        north = wrap(phase[1:] - phase[:-1])
        north[land[:-1] | land[1:]] = 0
        # swirl inside each little grid square, counted counter-clockwise (in full turns)
        swirl = (east[:-1] + np.roll(north, -1, 1) - east[1:] - north) / (2 * np.pi)

        touches_land = land[:-1] | land[1:] | np.roll(land[:-1], -1, 1) | np.roll(land[1:], -1, 1)
        water = np.rint(swirl[~touches_land]).astype(int)

        # group land cells into landmasses (diagonal neighbours and the date line join)
        lab, n = ndimage.label(land, structure=np.ones((3, 3)))
        parent = list(range(n + 1))

        def root(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(nlat):
            for k in (i - 1, i, i + 1):
                if 0 <= k < nlat and lab[i, -1] and lab[k, 0]:
                    parent[root(lab[i, -1])] = root(lab[k, 0])
        lab = np.vectorize(lambda v: root(v) if v else 0)(lab)

        coast = {}
        for k in np.unique(lab[lab > 0]):
            m = lab == k
            around = m[:-1] | m[1:] | np.roll(m[:-1], -1, 1) | np.roll(m[1:], -1, 1)
            coast[k] = int(round(swirl[around].sum()))

        print(f"\n{name.upper()}  (period {c['period_h']:.2f} h)")
        print(f"  still points in the water: {(water == 1).sum()} counter-clockwise, "
              f"{(water == -1).sum()} clockwise, net {water.sum():+d}")
        print(f"  swirl carried around coasts: net {sum(coast.values()):+d} "
              f"({sum(1 for v in coast.values() if v)} of {len(coast)} landmasses carry one)")
        print(f"  water + coasts = {water.sum() + sum(coast.values()):+d}   (topology says 0)")
        for place, (la, lo) in PLACES.items():
            i = int(round((la + 90) * 2))
            j = int(round(((lo + 180) % 360) * 2)) % nlon
            print(f"    {place:26s} {coast.get(lab[i, j], 0):+d}")


if __name__ == "__main__":
    main()
