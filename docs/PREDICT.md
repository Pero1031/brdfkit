# Predicting a plastic that was not measured, and checking it

hpBRDF has 14 materials, too few to learn a space of materials from. But its plastics are made of two parts that can be
taken apart by polarization ([UMOV.md](UMOV.md)): a surface reflection, white and polarized, and a body, coloured and
depolarized. If that is all there is, a plastic of a new colour needs no new table: a black plastic gives the surface, a
white one the body, and one spectrum of the new plastic says how much body there is at each wavelength. The data set
contains plastics that were not used to build this, so the prediction can be scored against a measurement.

```bash
python -m brdfkit predict out/quirks/Plum_rough_plastic.npz --white out/quirks/White_rough_plastic.npz \
       --black out/quirks/Black_rough_plastic.npz --glass out/quirks/Black_glass.npz
```

It works on the reductions of `brdfkit quirks` (mean Mueller matrices per θ_d and wavelength at the specular core,
θ_h ≤ 0.55°; beside the lobe, 15°–25°; far from it, ≥ 40°), so it runs in seconds. Implementation in
`brdfkit/analysis/predict.py`.

## The prediction

    M_T(geometry, λ) = M_B + r(λ) (M_W − M_B),      r(λ) = (R_T − R_B) / (R_W − R_B)

for all 16 Mueller elements: at every geometry and wavelength the target T lies between the black rough plastic B and the
white rough plastic W. The weight r comes from **one reading** of the target, its spectrum R far from the lobe at θ_d 30°,
which is what a colour measurement gives. The instrument's floor is in all three tables and drops out.

At the specular core the lobe has to be as wide as the target's, which the black plastic's (half width 8.3°) need not be.
There the white plastic's own surface is used (its m11: a surface reflection keeps the polarization, so m11 is its m00),
plus r times the white plastic's body under the peak. With **two more readings**, the half width w of the target's lobe and
its index n, that surface is scaled by (w_W / w_T)² F(n_T) / F(n_W): a lobe of the same energy and another width, and Fresnel.

The baseline is the **tint**, M_W × R_T / R_W: the white plastic multiplied by the colour, which is what colouring a
material in a renderer does. It dims the highlight with the colour and keeps the white plastic's polarization.

## Protocol

The method was fixed on the red rough plastic and committed (069128f) before the two targets, the yellow and the plum rough
plastic, were downloaded and reduced. The scores below are of that commit, unchanged. They cover θ_d 25°–75° in 5° steps
and 414–900 nm in the three regions (33 rows per target; the one the reading was taken from is among them).

## Result

"m00" is the median factor by which the prediction is off; "pol." the RMS error of m01, m10, m11, m22, m33 over m00.

| | | Red (development) | | Yellow (held out) | | Plum (held out) | |
|---|---|---|---|---|---|---|---|
| | Region | m00 | pol. | m00 | pol. | m00 | pol. |
| Tint | core | ×1.24 | 0.113 | ×2.01 | 0.138 | ×1.26 | 0.064 |
| | beside | ×1.41 | 0.145 | ×1.07 | 0.131 | ×1.25 | 0.196 |
| | far | ×1.23 | 0.029 | ×1.09 | 0.020 | ×1.09 | 0.065 |
| Interpolation, one reading | core | ×1.34 | 0.054 | ×1.69 | 0.095 | ×2.77 | 0.124 |
| | beside | ×1.23 | 0.042 | ×1.07 | 0.125 | ×1.07 | 0.038 |
| | far | ×1.08 | 0.023 | ×1.07 | 0.020 | ×1.09 | 0.047 |
| Interpolation, three readings | core | ×1.15 | 0.054 | ×1.42 | 0.163 | ×1.10 | 0.059 |

(Beside and far from the lobe the three-reading prediction is the one-reading one.)

- **The colour part works.** Away from the lobe, m00 of both held-out plastics is predicted within 7 %–9 % (median) over
  all angles and within 4 %–6 % up to θ_d 55°, from one spectrum. For the yellow plastic that is an extrapolation: its body
  is 1.5 times brighter than the white plastic's above 600 nm (r = 1.5), so the "white" plastic is not the brightest body
  there is, and the line through black and white still holds beyond white.
- **The plum plastic is predicted as a whole.** With the lobe's width and the index, its specular core is within 10 %, its
  polarization beside the lobe within 0.04 (the tint: 0.20), and the surface the model assumes is the one measured (m11
  beside the lobe 2.3 × 10⁻³ predicted, 2.4 measured at θ_d 30°).
- **The yellow plastic's surface is not the others'.** Its polarized part is 2.8 times the black plastic's beside the
  lobe and about twice what its lobe width and index predict at the core. Its colour is predicted, its highlight and its
  polarization are not (0.13 beside the lobe, no better than the tint; the three-reading core is worse than the
  one-reading one). The assumption that fails is "all these plastics have the same surface up to a width and an index".
- **Without the lobe's width the highlight is wrong** (plum ×2.8): the white plastic's lobe is 4.0° wide, the plum's 6.3°.
- Beyond θ_d 60° every prediction drifts (far from the lobe the plum is over by 1.2–1.4, the yellow by 1.2–1.6): the
  bodies do not share one dependence on θ_d there, which is the uneven cosine of [QUIRKS.md](QUIRKS.md) section 10.

### What it would take (after the fact)

Seen after the scores, so not a test: give the prediction two polarized readings of the target at θ_d 30°, its m11 at the
core and beside the lobe (a polarization-difference measurement), and scale the donors' surfaces by them
(`out/scripts/predict_posthoc.py`, not part of the package). Up to θ_d 55°:

| | Red | Yellow | Plum |
|---|---|---|---|
| Surface beside the lobe / black plastic's | 0.99 | 2.84 | 1.13 |
| m00: core, beside, far | ×1.08, ×1.20, ×1.05 | ×1.07, ×1.15, ×1.06 | ×1.04, ×1.09, ×1.04 |
| Polarization: core, beside, far | 0.050, 0.033, 0.021 | 0.047, 0.024, 0.012 | 0.057, 0.035, 0.014 |

With the surface's level known, all three are predicted to about the same accuracy. So the two-part structure holds on
every plastic checked; what does not transfer between samples is how strong the surface reflection is.

## What it means

- A coloured plastic's spectral, polarimetric reflectance is, to 5 %–20 % in m00 and 0.02–0.06 in the normalized Mueller
  elements, three things: a body shape shared with a white plastic, a spectrum, and a surface. The spectrum is one
  measurement. The surface is not predictable from another plastic; it needs its own measurement (width, level).
- Multiplying a material by a colour is wrong in a specific way that this data shows: it tints the highlight (the yellow
  plastic's core is off by 2) and it leaves the polarization that of the white material, where a dark colour polarizes
  strongly (plum beside the lobe: 0.20 against 0.04).
- For editing: recolouring a plastic should change r(λ) and leave the surface, which also changes the polarization the
  right way. `brdfkit swap` exchanges S and N as wholes, which is right for metals only.

## Limits

Three targets, all rough plastics from one data set, one white and one black donor. The scores are on three θ_h regions
summed over φ_d, not on single directions: the shape of the lobe between 0.55° and 15° is not tested. The readings are
taken from the target's own table (a spectrophotometer and a glossmeter would give them, but none was used). The lobe's
half width is measured on m00 above the body, which is not the surface's alone. Nothing is corrected for the quirks of the
tables; the comparison is between tables that share them.
