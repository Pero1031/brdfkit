# Where RGB stops being enough: inter-reflection between the measured plastics

A renderer that carries three numbers per colour multiplies the light's RGB by a surface's RGB at every bounce. The spectrum
is multiplied per wavelength and turned into a colour at the end. The two agree while the spectra are smooth. The dyed
plastics of hpBRDF have steep absorption edges, and this page measures how far the RGB product drifts on them: after one
reflection under daylight it is right (ΔE00 about 1), after two it is off by 6 to 7 on the red and the yellow plastic,
under a tungsten lamp or a screen's light it is off from the first reflection, and a wider gamut does not repair it.

```bash
python -m brdfkit bounce out/quirks/Red_rough_plastic.npz out/quirks/Yellow_rough_plastic.npz out/quirks/Plum_rough_plastic.npz \
       out/quirks/White_rough_plastic.npz --white out/quirks/White_smooth_plastic.npz --reference out/quirks/Black_glass.npz
```

It works on the reductions of `brdfkit quirks`. A sample's reflectance ρ(λ) is its intensity beside the lobe (θ_h 15°–25°,
θ_d 30° and 45°, floor removed) relative to the white smooth plastic's, times an albedo given to that white (0.65),
414–774 nm. A surface's RGB is the colour of ρ under an equal-energy light, a light's RGB the colour of its spectrum; the
RGB result after n reflections is light × surfaceⁿ per channel. Both results are brought to the same exposure (the spectral
one to a mid grey) before the difference is taken, so ΔE00 is what a viewer of the two pictures would see, and the number in
brackets is the luminance of the RGB result over the spectral one. Lights: a 6500 K and a 2856 K blackbody, and a screen
with three bands (455, 535, 625 nm) balanced to the daylight's colour. Implementation in `brdfkit/analysis/bounce.py`.

## The same plastic, n reflections

| Linear sRGB, daylight | n = 1 | 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|
| Red rough | 1.0 (0.97) | 5.6 (1.13) | 6.5 (1.24) | 6.8 (1.29) | 6.7 (1.29) |
| Yellow rough | 1.3 (0.99) | 7.4 (1.00) | 13.3 (1.05) | 16.0 (1.16) | 19.6 (1.59) |
| Plum rough | 0.4 (0.99) | 0.3 (1.00) | 0.5 (1.02) | 1.8 (1.06) | 5.4 (1.17) |
| White rough | 0.0 (1.00) | 0.1 (1.00) | 0.3 (1.00) | 0.6 (1.00) | 1.1 (1.01) |

One reflection is what RGB was made for, and it holds. The second reflection off the red or the yellow plastic is already
a visible error. The plum plastic, whose spectrum is nearly flat, holds for three reflections, the white for any number.
For comparison, gold needs three reflections to reach 3.7 (GoLD_rayt, `docs/METAL_ACCURACY.md`): a dye's edge is steeper
than a metal's.

The reason is in the edge. A channel's value is the reflectance averaged over the channel; the red plastic's red channel
averages over its edge at 600 nm. After several reflections the light that is left sits where the reflectance is highest,
beyond the edge, and that average no longer describes it.

## Other lights

| Linear sRGB | n = 1 | 2 | 3 |
|---|---|---|---|
| Red, incandescent | 8.5 (1.09) | 15.0 (1.26) | 15.1 (1.35) |
| Yellow, incandescent | 8.7 (1.04) | 18.4 (1.11) | 22.8 (1.23) |
| Plum, incandescent | 2.0 (1.03) | 4.2 (1.07) | 6.8 (1.12) |
| Red, three-band screen | 7.3 (0.91) | 6.3 (0.87) | 7.1 (0.84) |
| Yellow, three-band screen | 5.4 (1.08) | 4.8 (1.19) | 7.8 (1.33) |
| Plum, three-band screen | 6.3 (1.03) | 7.5 (1.05) | 7.4 (1.07) |
| White, either | ≤ 0.4 | ≤ 0.8 | ≤ 1.1 |

Under a light that is not the one the RGB values were made for, the first reflection is already wrong by 5 to 9 on the
dyed plastics: the product of two averages is not the average of the product. The screen's light has the same colour as the
daylight and gives a different result, because its three bands sample the spectrum at three places only; a room lit by a
monitor is a case for a spectral renderer, the more so with the screen's polarization
(GoLD_rayt, `docs/POLARIZATION.md`).

## A wider gamut does not repair it

| Linear Rec. 2020, daylight | n = 1 | 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|
| Red rough | 0.9 (1.03) | 7.0 (0.75) | 12.8 (0.54) | 17.7 (0.39) | 26.4 (0.19) |
| Yellow rough | 0.5 (1.01) | 3.5 (0.94) | 8.1 (0.86) | 12.0 (0.78) | 16.5 (0.67) |
| Plum rough | 0.3 (1.01) | 0.7 (0.98) | 1.7 (0.93) | 3.2 (0.86) | 7.5 (0.70) |

Narrower primaries help the yellow plastic and hurt the red one, and they turn the error from too bright into too dark:
after six reflections the red plastic keeps 0.19 of the right luminance. The error is not one of gamut; three numbers
cannot hold where inside a channel the reflectance sits.

## Two different plastics

One reflection off each, linear sRGB, daylight: red then yellow 6.8 (1.10); the other pairs 0.3 to 1.9. Under the
incandescent light: 8.1 to 16.5 for every pair with the red or the yellow plastic, 1.8 for plum then white.

## Check with a renderer

An open box whose five walls are one plastic (a diffuse surface with ρ(λ)), a small daylight panel at its mouth, rendered
by the spectral and the RGB build of GoLD_rayt from the same scene file (`out/scripts/bounce_scene.py`, 480 px, 256 spp).
ΔE00 between the two pictures at equal exposure, and the luminance of the RGB picture over the spectral one:

| Box | whole box | darkest third (mostly bounced light) |
|---|---|---|
| White rough | 0.2 (1.00) | 0.3 (1.01) |
| Red rough | 1.3 (1.06) | 3.1 (1.15) |
| Yellow rough, ρ up to 0.97 | 21.2 (1.87) | 25.0 (2.21) |
| Yellow rough, ρ scaled to at most 0.75 | 10.2 (1.18) | 14.2 (1.32) |

The yellow box comes out a redder orange and brighter in RGB. The red box differs less than the table of reflections
suggests, because the red plastic is dark (ρ at most 0.4): most of what the camera sees has bounced once. The yellow box
depends on the level given to its reflectance, which these tables do not fix (below); between the two levels tried the
error is 10 to 21.

## Limits

- The reflectance is relative: the tables give each sample's intensity against the white sample's, and the white's albedo
  is known only as 0.68–0.75 ([QUIRKS.md](QUIRKS.md)). With 0.7 the yellow plastic would reflect more than all the light at
  600 nm, so the level of one sample against another is uncertain by tens of percent. The tables of reflections do not
  depend on it (each result is exposed separately); the box does.
- ρ is read beside the lobe, where a rough sample's surface reflection adds a little white light; the body alone has
  steeper edges.
- The surfaces are taken as diffuse with one spectrum. The surface reflection, which is white and would dilute the colour
  at each bounce, is left out; so is the fluorescence some yellow dyes have.
- The colour matching functions are an analytic fit, the daylight a blackbody, the screen three Gaussians.
- An RGB renderer can be made better than the plain product (spectral upsampling of each RGB); the comparison here is with
  the plain product, which is what most do.
