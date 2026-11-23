# Other eyes: a near-infrared camera and a polarization camera on the measured plastics

The hpBRDF tables hold 68 wavelengths from 414 to 950 nm and the full Mueller matrix, which is more than an eye uses. This
page asks what two other eyes would see of the six plastics, away from the highlight and under unpolarized light: a camera
with a near-infrared filter, and a polarization camera. The answers are that the near infrared removes the colours of the
dyed plastics but not all of them, and that the polarization a camera sees on a plastic depends on the band so strongly that
the choice of filter decides whether there is anything to see: on the red plastic 60 % in green and 0.4 % in red.

```bash
python -m brdfkit eyes out/quirks/Black_rough_plastic.npz out/quirks/Plum_rough_plastic.npz out/quirks/Red_rough_plastic.npz \
       out/quirks/Yellow_rough_plastic.npz out/quirks/White_rough_plastic.npz out/quirks/White_smooth_plastic.npz \
       --reference out/quirks/Black_glass.npz --white out/quirks/White_smooth_plastic.npz
```

It works on the reductions of `brdfkit quirks`: Mueller matrices summed over φ_d at θ_h 15°–25°, beside the lobe. The
instrument's unpolarized floor is taken from the black glass in the same region ([FLOOR.md](FLOOR.md)). Bands: blue
430–490 nm, green 510–570, red 600–660, NIR 800–900. Implementation in `brdfkit/analysis/eyes.py`.

## What a camera records

Under unpolarized light a camera gets the first column of the Mueller matrix: the intensity R = m00 and the linearly
polarized part q = m10 (m20 is small here). Two kinds of light reach it ([UMOV.md](UMOV.md)):

- light reflected by the surface: not coloured by the pigment, polarized along s (q > 0 in the table's frames);
- light that went through the body: the pigment's colour, and weakly polarized along p (q < 0) by the refraction on its
  way out.

So q(λ) = q_surface − |q_body| B(λ), with B the body's spectrum. The surface's share is the same at every wavelength; the
body's follows the colour.

## The near infrared: which colours go

Intensity relative to the white smooth plastic, θ_d 45°:

| | blue | green | red | NIR |
|---|---|---|---|---|
| Black rough | 0.03 | 0.05 | 0.05 | 0.09 |
| Plum rough | 0.11 | 0.13 | 0.22 | 0.24 |
| Red rough | 0.04 | 0.05 | 0.63 | 0.70 |
| Yellow rough | 0.23 | 0.78 | 1.46 | 1.42 |
| White rough | 1.10 | 1.03 | 1.01 | 1.02 |

The red and the yellow plastic are as light in the near infrared as in the red band: their pigments stop absorbing at
600 nm and 540 nm and do not start again up to 950 nm. Through an infrared filter the red plastic is a light grey, 0.70 of
the white, where the eye's green sees 0.05; red on white loses most of its contrast. The plum plastic does not follow: it
stays at 0.24, flat from 600 to 950 nm, so its darkness is not that of a red dye but of something that absorbs at all
wavelengths. The black stays black. An infrared camera therefore sorts these samples by what darkens them, not by hue.
(The yellow is brighter than the white sample, whose albedo is only about 0.7.)

## Polarization: the band decides

q / R in percent (+ along s, the surface; − along p, the body), with q × 1000 in brackets; "dark" marks an intensity below a
third of the floor, where the ratio is unreliable:

| θ_d 45° | blue | green | red | NIR |
|---|---|---|---|---|
| Black rough | dark (2.5) | 89 (2.4) | 82 (2.5) | 65 (3.0) |
| Plum rough | 26 (1.7) | 24 (1.8) | 15 (1.8) | 16 (2.0) |
| Red rough | dark (1.8) | 60 (1.8) | 0.4 (0.1) | 0.5 (0.2) |
| Yellow rough | 28 (3.7) | 5.9 (2.6) | 0.1 (0.1) | 1.4 (1.0) |
| White rough | −3.1 (−2.0) | −2.4 (−1.4) | −1.9 (−1.1) | 0.7 (0.4) |
| White smooth | −6.6 (−3.8) | −6.1 (−3.5) | −5.7 (−3.2) | −3.3 (−1.7) |

- **Where the pigment absorbs, the surface is all there is**, and the light is strongly polarized: the red plastic 60 % in
  green. Where the pigment lets go, the body adds a hundred times more unpolarized light and its own opposite
  polarization: 0.4 % in red and in the near infrared. A polarization camera with a red or infrared filter sees nothing on
  it; with a green or blue one it sees the surface.
- **The polarized part itself, q, is the same in every band on the dark samples** (black 2.4–3.0, plum 1.7–2.0): it is the
  surface, and does not know the colour. On the red and the yellow plastic it falls from 1.8 and 3.7 to 0.1 because the
  body's negative share cancels it.
- **On the white plastics the sign is the body's.** The angle of polarization is turned by 90° against the coloured
  samples in their dark bands. Methods that read the surface's orientation from the angle of polarization have to know
  which of the two it is; on a coloured plastic the two bands give both.
- **The sign changes inside the spectrum** where the two shares are equal. The yellow plastic at θ_d 30° goes from
  +6.6 % in blue to −0.9 % in red, through zero at 562 nm, at its absorption edge.

So the rule for a polarization camera on a dyed plastic is the opposite of the usual one for a dark scene: use the band in
which the object is darkest. The near infrared is the worst band for the dyed samples and no worse than any other for the
black and the plum one.

## Polarization: glossy or matt, without the highlight

The two white plastics have the same colour and, away from the highlight, the same intensity (within 10 %). Their
polarization differs, and more so at larger angles:

| q / R [%], green | θ_d 30° | 45° | 60° | 70° |
|---|---|---|---|---|
| White smooth | −2.3 | −6.1 | −11.1 | −6.9 |
| White rough | −1.6 | −2.4 | +2.4 | +22.2 |

On the smooth sample the surface reflects into the highlight only, so beside it the camera sees the body alone (p). On the
rough one the surface's lobe is wide and reaches these directions, and at 60° and beyond it outweighs the body (s). The
sign of the polarization, 15° to 25° away from the mirror direction, tells a glossy white from a matt white without the
highlight being in the picture.

## No single band

The weakest |q| / R among the six samples, per band:

| θ_d | blue | green | red | NIR |
|---|---|---|---|---|
| 30° | 2.3 | 1.1 | 0.9 | 0.2 |
| 45° | 3.1 | 2.4 | 0.1 | 0.5 |
| 60° | 1.3 | 2.4 | 3.4 | 6.1 |
| 70° | 4.8 | 6.9 | 2.7 | 2.9 |

No band keeps every sample above a few percent: the whites are weak everywhere, and each dyed sample has bands where it
vanishes. Up to 45° blue is the safest, because every pigment here absorbs it or, on the whites, the body's share is
largest there.

## Limits

- Six samples. The statements are about these pigments; the near-infrared behaviour of the plum one is read from its
  spectrum, not from knowing what is in it.
- 950 nm is the end of the data. Plastics are sorted by their absorption bands beyond 1000 nm, which these tables do not
  reach.
- The geometry is one region beside the lobe, summed over φ_d, in the table's frames; the sums weigh directions by the
  cosine of the light's angle that the tables carry ([QUIRKS.md](QUIRKS.md), 10). Ratios within one sample are not
  affected, comparisons between samples only through their different lobes.
- The floor is the black glass's, with scale 1. On the black plastic the intensity left after removing it is of the size of
  its uncertainty (ratios above 100 % occur at θ_d 60°); its q is reliable, its q / R is not.
- The whites' polarization fades toward the near infrared (−6.6 % to −3.3 % on the smooth one) while the black plastic's q
  is flat. This is as measured and not explained here.
- Metals and the glass are left out: beside the lobe they return the floor and nothing else. Their specular polarization
  is in [QUIRKS.md](QUIRKS.md) and [ELLIPSO.md](ELLIPSO.md).
