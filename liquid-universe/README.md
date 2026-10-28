# Amphidromic Points Explained Visually

Two interactive models that show what **amphidromic points** are, using a real map of the ocean tide.

![Where the Tide Stands Still](youtube_thumbnail_barrel.png)

## What is an amphidromic point?

It's a spot in the ocean where the sea level barely rises or falls, even though the tide is moving all around it.

The tide doesn't just slosh in and out. In each ocean basin it travels around in a big circle, so high tide sweeps around one side while low tide sweeps around the other. Right in the middle, the two balance out and the water stays at about the same level. It's like the hands of a clock: they keep turning, but the centre they turn around never moves.

There are dozens of these points around the world. Most of them turn counter-clockwise in the Northern Hemisphere and clockwise in the Southern, because of Earth's rotation.

## Try the models (free, runs in your browser)

| Model | What it shows | Live link |
|---|---|---|
| **Globe** | The real ocean tide (red = high, blue = low) with every major amphidromic point marked. Click anywhere to see a cross-section of the water, a tide gauge, and how the current underneath flips when the tide turns. | [amphidromic_earth.html](https://krimseylilleth.github.io/universe/liquid-universe/amphidromic_earth.html) |
| **Close-up** | A simple "kiddie pool" with one flat wave going around it, showing how a still point in the middle is possible. Also a wide-area view of a tide wave travelling around a channel. | [amphidromic_closeup.html](https://krimseylilleth.github.io/universe/liquid-universe/amphidromic_closeup.html) |

## Video

*Amphidromic Points Explained Visually*: link coming soon.

## Controls

- Drag to rotate, scroll to zoom
- Click the ocean (globe) to probe the tide at that spot
- **Press S** to hide / show all panels
- The globe can switch between tide constituents: **M2** (the Moon's twice-a-day tide, 12.42 h), **S2** (the Sun's, 12.00 h), and the once-a-day **K1** and **O1**

## Data

- **Tides:** NASA GOT4.10c ocean tide model (R. Ray, NASA Goddard), built from satellite altimetry since 1992, on a 0.5° grid. <https://earth.gsfc.nasa.gov/geo/data/ocean-tide-models>
- **Sea floor depth:** NOAA ETOPO1, averaged to 0.5°. <https://www.ncei.noaa.gov/products/etopo-global-relief-model>
- Currents are calculated from the slope of the tide surface and Earth's rotation (linear, frictionless), not measured. See the "Honest notes" in the globe's panel for the limits.

## Files

- `amphidromic_earth.html`: the globe
- `amphidromic_closeup.html`: the close-up (self-contained)
- `amphidromic_got410c_data.js`, `amphidromic_bathy_data.js`: tide and depth data used by the globe

## Suggestions welcome

Ideas or fixes? [Open an issue](https://github.com/KrimseyLilleth/universe/issues), or click the pencil icon on any file to propose an edit.

Thumbnail photo: Austin Schmid via Unsplash (CC0).

Krimsey Lilleth · klilleth@proton.me
