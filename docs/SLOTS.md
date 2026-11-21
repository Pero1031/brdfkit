# Which direction is the light? The hpBRDF tables are filed at φ_d + 180°

A table of Mueller matrices over (φ_d, θ_d, θ_h) names two directions, and a reader has to know which one is the light. The
hpBRDF tables are read with the Mitsuba `measured_polarized` plugin, which takes the first Rusinkiewicz direction as the light.
The polarization of a diffuse body shows that in the hpBRDF tables the first direction is the camera: the table holds the
right Mueller matrices, filed at φ_d + 180°. This page gives the check, its result, and the correction.

```bash
python -m brdfkit slots data/kaist/hpbrdf/White_smooth_plastic.hpbrdf --out out/slots/White_smooth_plastic.npz
python -m brdfkit slots data/kaist/hpbrdf/Red_rough_plastic.hpbrdf --bands 620,720 --theta-h-min 35 --out out/slots/Red_rough_plastic.npz
python -m brdfkit slots out/slots/White_smooth_plastic.npz --out out/slots/White_smooth_plastic.npz     # report again, in a second
```

It reads every sixth φ_d slice (about two minutes for a 13 GB file). Implementation in `brdfkit/analysis/slots.py`.
"As the plugin labels them" is not an assumption of this code: tables written by brdfkit were read back with the stock
plugin of Mitsuba 3.7.1 and agree with brdfkit's own reading in all 16 Mueller elements ([MITSUBA.md](MITSUBA.md)), and the
authors' plugin in the hpBRDF repository calls `directions_to_rusinkiewicz` with the same order of arguments.

## The check

Light scattered inside a dielectric leaves through its surface. Whatever went in, what comes out is partly polarized by the
transmission through the interface: in the plane of the exit direction and the normal (p), and more so the more grazing the
exit. In the same way, how much light gets in depends on its polarization through the angle it arrives at. In a Mueller
matrix (S_out = M S_in) the first is the first column, (m10, m20), and the second the first row, (m01, m02). So on a
diffuse body:

- |(m10, m20)| / m00 is a function of the **camera's** angle alone, and the vector points along the camera's plane of emergence;
- |(m01, m02)| / m00 is a function of the **light's** angle alone, and points along the light's plane of incidence.

This needs no reference sample and no calibration, only a sample with a body under a smooth interface. The command bins the
body (θ_h ≥ 15°, away from the specular lobe; θ_d 22°–78°, the measured range) by the two polar angles as the plugin labels
them, and reports

1. how much of the variance of each magnitude one angle alone explains;
2. the measured vectors projected on the predicted orientation, for three readings of the table: as labelled; the two
   directions exchanged; exchanged with the signs of S2 and S3 flipped ("mirrored", a frame of the other handedness);
3. the exponents of m00 = c cos^a("light") cos^b("camera").

## Result

White smooth plastic, 550–694 nm. Polarizance |(m10, m20)| / m00 by the plugin's two angles (rows: "light", columns: "camera"):

| "light" \ "camera" | 0°–10° | 20°–30° | 40°–50° | 60°–70° | 70°–80° |
|---|---|---|---|---|---|
| 0°–10° | – | – | 0.002 | 0.005 | 0.008 |
| 20°–30° | – | – | 0.012 | 0.013 | 0.015 |
| 40°–50° | 0.044 | 0.047 | 0.048 | 0.045 | 0.039 |
| 60°–70° | 0.125 | 0.125 | 0.127 | 0.123 | 0.103 |
| 70°–80° | 0.186 | 0.187 | 0.182 | 0.171 | 0.132 |

The polarization of the light that comes out follows the angle the plugin calls the light's, and not the one it calls the
camera's: the "light" angle alone explains 0.97 of the variance, the "camera" angle 0.07. The values are those of
transmission out of a dielectric (index 1.62 from this sample's specular polarization: 0.057 at 45°, 0.24 at 75°). The
sensitivity to the incoming polarization does the opposite (0.01 against 0.93).

| Reading | Column (m10, m20) on the prediction | its S2 term | Row (m01, m02) | its S2 term |
|---|---|---|---|---|
| as labelled | −0.13 | −0.70 | −0.36 | −0.79 |
| **exchanged** | **+0.99** | **+1.00** | **+0.79** | **+0.99** |
| exchanged, mirrored | −0.03 | −1.00 | −0.18 | −0.99 |

Only one reading agrees with the orientation as well: the directions exchanged, with each beam keeping the plugin's
reference vector and no change of sign. The S2 term, which is zero in the plane of incidence and so tests the frames out of
the plane, agrees to 1.00.

The red rough plastic (620–720 nm, where its body is bright) gives the same reading once the rough surface's reflection is
left behind at θ_h ≥ 35°: polarizance explained by the "light" angle 0.97; orientation +0.98 (S2 term +1.00) exchanged,
−0.22 as labelled, −0.26 mirrored. At θ_h ≥ 15° the surface reflection, which is s-polarized, still cancels part of it
(+0.46, S2 term +0.99).

## What it means

The plugin's reference vectors only change sign when the two directions trade places (x_in(a, b) = −x_out(b, a), checked in
the tests), and a sign of x does not change a Stokes vector. So "exchanged, not mirrored" means: the matrix stored at
(φ_d, θ_d, θ_h) is the correct matrix, in the plugin's own frames, for the geometry with light and camera exchanged, which
is the bin at φ_d + 180°. The correction is a shift of the φ_d axis by 180° and nothing else:

```python
from brdfkit.analysis import clean
table = clean.exchange(table)          # or --exchange in compress, eval, factor, slots, to-mitsuba
```

`brdfkit slots TABLE --exchange` on the corrected view reads "as labelled".

What changes and what does not:

- A specular reflection is the same with light and camera exchanged, so everything measured on the specular lobe at θ_h ≈ 0
  (lobe energy, calibration, the polarization against Fresnel in [QUIRKS.md](QUIRKS.md)) is unaffected.
- The polarization of a diffuse body was attached to the wrong direction: rendered as stored, the polarization pattern of a
  plastic follows the light's angle where it should follow the camera's. To the unpolarized image it makes no difference.
- The one-sided cosine in m00 (quirk 10) is on the side the plugin calls the camera's. With the directions exchanged that is
  the light: the side on which the paper's fit, which has no cosine of the irradiance, puts it. One exchange accounts for
  the polarization and for the side of the cosine.

## The cosine is not the same on every sample

The exponents of m00 = c cos^a("light") cos^b("camera") on the body, as stored (no floor removed, all angles below 80°):

| | a ("light": the camera) | b ("camera": the light) |
|---|---|---|
| White smooth plastic, θ_h ≥ 15° | −0.04 | +0.78 |
| Red rough plastic, θ_h ≥ 15° | −0.15 | +0.23 |
| Red rough plastic, θ_h ≥ 35° | −0.08 | +0.27 |

A reciprocal body has a = b, so b − a is what the table adds. The white plastic carries most of one cosine of the light
(more with the floor removed and the angles kept below 45°: 1.1–1.2, [QUIRKS.md](QUIRKS.md) section 10); the red rough
plastic carries only about a third of one. So dividing by the light's cosine (`--cosine light`, after `--exchange`) is right
for the white plastics and for the lobes of the smooth samples, whose reciprocity it restores, and too much for the red one.
The tables do not hold one factor for all samples; `slots` prints the exponents of a sample so the choice can be made per
table.

## Limits

- Decided on two samples, the white smooth and the red rough plastic; both give the same reading. The white rough plastic's
  table was not kept. Metals and glass have no body and cannot be checked this way; their specular data do not depend on it.
- The check assumes that a body polarizes like light leaving a dielectric. A table that does not is left undecided: the
  KAIST pBRDF 2020 Spectralon target (sintered PTFE, no interface) gives 0.29 against 0.37 and no reading above 0.5.
- The sign of S3 (circular polarization) is not tested by a body, which produces none. A frame of the other handedness
  flips S2 and S3 together and S2 is settled, so S3 is taken to be right as stored.
- A table stored as the transpose of its matrices at the labelled geometry would look the same on a reciprocal sample.
  The side of the cosine in m00 speaks for the exchange, since the paper's fit predicts it on the light's side.
