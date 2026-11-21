# Surface and body, split by polarization

A pigmented plastic reflects in two ways: off its surface, and from inside. The two differ in colour and in polarization,
and hpBRDF measures both, so on a coloured sample they can be taken apart without a model of either. This page does it for
the red rough plastic, with the black rough plastic as a control. It is a check of the data (the dielectric samples'
polarization and their spectra away from the specular lobe behave as physics says) and a tool: the surface reflection of a
plastic, lobe and tail, read out from under its body.

```bash
python -m brdfkit quirks data/kaist/hpbrdf/Red_rough_plastic.hpbrdf --out out/quirks/Red_rough_plastic.npz --dielectric
python -m brdfkit quirks data/kaist/hpbrdf/Black_glass.hpbrdf --out out/quirks/Black_glass.npz --dielectric
python -m brdfkit umov out/quirks/Red_rough_plastic.npz --reference out/quirks/Black_glass.npz --n 1.505
```

`umov` works on the reductions that `brdfkit quirks` writes (Mueller matrices summed over φ_d, beside the lobe at
θ_h 15°–25°, `--region mid`, or far from it at θ_h ≥ 40°, `--region far`), so it runs in a second. Implementation in
`brdfkit/analysis/umov.py`.

## The model

At one geometry, over wavelength:

    M(λ) = M_surf + B(λ) N_body + k c(λ) e00

- **M_surf**: Fresnel reflection off the facets of the rough surface. It is white (the index of a plastic hardly changes
  over the range) and keeps the polarization: in the table's frames m11 = m00.
- **B(λ) N_body**: light scattered inside. B is the pigment's spectrum; N_body is close to e00, a depolarizer.
- **k c(λ) e00**: the instrument's unpolarized floor ([FLOOR.md](FLOOR.md)). Its spectrum c is taken from black glass in the
  same region (m00 − m11 of the glass), and its scale k is fitted on the bands where the pigment absorbs everything
  (below 540 nm for the red plastic): there m00 − m11 = k c + u, with u a white unpolarized remainder.

With R = m00 − k c the reflection without the floor, the share of polarized light is m11 / R, and since m11 belongs to the
surface and does not change with wavelength, it falls as 1 / R. That is Umov's law (1905): the darker, the more polarized.
Here it is tested on one sample, over the 14-fold change of R that the pigment makes between 550 and 620 nm.

## Result on the red rough plastic (θ_h 15°–25°)

| θ_d | floor k | surface m11 × 10³ | m01 / m11 (Fresnel) | m22 / m11 (Fresnel) | body m11, m10 | R max / min | Umov slope |
|---|---|---|---|---|---|---|---|
| 30° | 0.84 | 2.17 | 0.41 (0.39) | −0.96 (−0.92) | +0.01, 0.00 | 14.3 | −0.93 |
| 40° | 0.84 | 2.23 | 0.62 (0.69) | −0.78 (−0.73) | 0.00, −0.04 | 14.5 | −1.00 |
| 45° | 0.91 | 2.41 | 0.78 (0.83) | −0.57 (−0.56) | −0.02, −0.05 | 14.3 | −1.09 |
| 50° | 0.90 | 2.82 | 0.90 (0.94) | −0.35 (−0.34) | +0.01, −0.06 | 13.6 | −0.97 |
| 55° | 0.85 | 3.84 | 0.92 (1.00) | −0.10 (−0.08) | +0.01, −0.08 | 10.7 | −0.96 |
| 60° | 0.89 | 5.45 | 0.94 (0.98) | +0.20 (+0.20) | +0.04, −0.09 | 7.9 | −0.89 |
| 70° | 0.90 | 16.3 | 0.75 (0.75) | +0.69 (+0.66) | +0.17, −0.01 | 3.8 | −0.73 |

Fresnel is that of the index fitted from the sample's own specular polarization (1.505), at the facet's angle of incidence,
θ_d. "Umov slope" is the slope of log(m11 / R) against log R over wavelength (−1 is the law).

- **The polarized part does not follow the colour.** R changes 14-fold over wavelength; m11, m22 and m33 stay.
- **It is a Fresnel reflection**, 15°–25° away from the specular direction: at θ_d 40°–70° m22 / m11 follows Fresnel within
  0.05 through the change of sign near Brewster's angle, m01 / m11 reads up to 0.08 low, and m01 ≈ m10 within 0.04 (m10
  reads 0.1 low below 35°, as it does at the specular of the glass, [CALIBRATION.md](CALIBRATION.md)).
- **Umov's law holds** at θ_d 25°–60°: slopes −0.85 to −1.09. From 65° the surface term grows, R varies by less than 6×
  and the slope falls off (−0.80 to −0.49).
- **The body depolarizes almost fully**: m11, m01, m22 and m33 of N_body are at or below 0.04 at θ_d 25°–55°. The exception
  is m10, −0.04 to −0.09 at 40°–65°: the sign opposite to the surface reflection, as for light that leaves through an
  interface ([SLOTS.md](SLOTS.md) uses this on the white plastic).
- Under unpolarized light the degree of polarization m10 / R at θ_d 50° is 0.78 at 502 nm and 0.007 at 702 nm; at 60°,
  0.88 and 0.05.

The black rough plastic, which has no body, is the control: its surface term is nearly the red one's (m11 × 10³ = 2.05
against 2.17 at θ_d 30°) and follows Fresnel of its own index (1.457) the same way (m22 / m11 within 0.04 at 45°–70°). The red plastic's surface, read out from
under a body 14 times brighter, is the surface of a plastic with nothing under it.

## What it is for

- **A check of the tables.** Away from the specular lobe, the dielectrics' polarization and spectra are consistent with a
  two-component reflection to a few percent. The quirks of [QUIRKS.md](QUIRKS.md) are in the scale and in the angles, not here.
- **The surface reflection of a diffuse material.** m11 is the surface, m00 − m11 − floor the body. For a rough plastic the
  lobe's energy cannot be measured on m00 (quirk 5: the lobe reaches into the body), but it can on m11.
- **Editing.** For a plastic the colour and the polarization are not independent: a darker colour is a more polarized one.
  `brdfkit swap` moves S and N separately, which is right for metals and for the specular lobe, not for a body.
- **The floor.** k comes out at 0.80–0.91 of the glass's floor on the red plastic and 0.57–0.86 on the black rough one:
  the floor differs from sample to sample by more than the ±15 % read from the far field alone.

## Limits

One coloured sample (the red rough plastic); the yellow and plum rough plastics of the data set would repeat it. The floor
fit assumes the body is black below `--dark-max`; the polarized elements do not depend on the floor at all, so the slopes
move by ±0.04 when k is changed by ±20 %. The regions are sums over φ_d, so the split is per (θ_d, θ_h range), not per
direction. Far from the lobe (`--region far`) the surface term is below 10⁻³ and the body's own residual polarization, a few
percent of a much brighter term, is as large: the ratios to Fresnel do not hold there. A white unpolarized remainder u of
0.1–0.4 × m11 is left in the dark bands (multiple reflection between facets, or body light that the pigment does not fully
absorb; not separated). The law itself is old; what is new is only that this data set shows it on one sample.
