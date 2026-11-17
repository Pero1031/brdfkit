# Gloss: one physical knob for the width of the lobe

Implementation in `brdfkit/analysis/gloss.py`.

```bash
# width of a lobe (hwhm of D at theta_d 20 / 30 / 50 / 70 deg, and its GGX-equivalent roughness)
python -m brdfkit gloss out/swap/Aluminum_clean.npz
# a rougher aluminium: to a GGX-equivalent roughness, or to a half width in degrees
python -m brdfkit gloss out/swap/Aluminum_clean.npz --alpha 0.21 --out out/gloss/Aluminum_a0.21.npz
python -m brdfkit gloss out/swap/Aluminum_clean.npz --hwhm 5 --out out/gloss/Aluminum_h5.npz
```

The input is a factor model, ideally the clean one of [SWAP.md](SWAP.md) (calibrated and symmetrized over φ_d). The result
is an ordinary factor model: `compress` and `export` bake it for GoLD_rayt, and `swap` can change its optics.

## What the knob does

A microfacet reflection is f = F(θ_d) D(θ_h) G(i, o) / (4 cos θ_i cos θ_o). D, the distribution of microfacet normals,
holds the roughness; F, the reflectance, and the polarization depend on θ_d only. In a factor model
([FACTOR.md](FACTOR.md)) every row of A at fixed (φ_d, θ_d) is therefore, up to the constant F(θ_d),

    D(θ_h) = A · 4 cos θ_i cos θ_o        (G ≈ 1 for the polished hpBRDF materials)

and P(m) = D cos⁴θ_h is the density of the microfacet slopes m = tan θ_h over the slope plane. A rougher surface of the
same material has the slopes of the original plus independent slope noise, so its slope density is the 2-D convolution of P
with a kernel: the slope density of GGX (default) or Beckmann of some roughness. Then

    A' = D' G' / (4 cos θ_i cos θ_o) = D' / (k_i k_o),   k = cos θ + sqrt(cos²θ + α'² sin²θ)

with G' the Smith shadowing of GGX at the result's alpha_eq α'. G' vanishes like cos θ_i cos θ_o toward the horizon, so A'
stays finite there; without it, the near-horizon bins of a widened lobe become arbitrarily bright (fireflies in a render).

- F, S and N are untouched, so are the colour and the polarization.
- The convolution keeps the energy of the lobe. Light that the wider lobe sends below the horizon, or that the shadowing
  removes, is lost, as for a real rough surface in single scattering. `--keep-albedo DEG` rescales A to keep the albedo at
  one angle of incidence instead.
- What the convolution only smooths stays: the haze tail of the material, and the instrument's response (next section).
- The knob is the width of the result: the half width at half maximum (hwhm) of the φ_d-averaged D(θ_h) at θ_d = 30°
  (`--ref-theta-d`), or alpha_eq, the GGX roughness whose D has that hwhm (`--alpha`). The kernel that gives it is found by
  bisection.

Checked on synthetic lobes (`tests/test_brdfkit.py::test_gloss`): Beckmann convolved with Beckmann gives Beckmann with the
roughnesses added in quadrature within 1% of the peak, the slope density keeps its mass within 2%, the result is bounded by
D'/α'² at the horizon, and a convolved GGX lobe has the directional albedo of the analytic GGX lobe with Smith shadowing
within about 1%.

## The measured lobes: instrument response convolved with the surface

Width of the four hpBRDF materials (calibrated, symmetrized):

| Material | hwhm at θ_d 20° / 30° / 50° / 70° | alpha_eq | Kernel that turns black glass into it |
|---|---|---|---|
| Black glass | 1.09° / 1.15° / 1.49° / 2.46° | 0.031 | — |
| Fake gold | 1.11° / 1.14° / 1.44° / 2.31° | 0.031 | none (as narrow as the glass) |
| SUJ2 | 1.45° / 1.59° / 1.88° / 2.80° | 0.043 | GGX 0.019 |
| Aluminum | 1.74° / 1.80° / 2.10° / 3.45° | 0.049 | GGX 0.027 |

Black glass is optically flat, so its lobe is the instrument's response: about 1.1° at θ_d 20°–30°, widening to 2.5° at 70°,
and wider in the plane of incidence than across it. Every material widens with θ_d in the same way, which a surface
cannot do (D depends on θ_h only). Fake gold is as narrow as the glass: within the instrument's resolution it is a mirror.

If a measured lobe is the glass lobe convolved with a roughness kernel, one kernel fitted at θ_d 30° should give the width at
every θ_d. It does, up to 50°–60°:

| hwhm at θ_d | 20° | 30° | 40° | 50° | 60° | 70° | 75° |
|---|---|---|---|---|---|---|---|
| SUJ2, measured | 1.45 | 1.59 | 1.70 | 1.88 | 2.15 | 2.80 | 3.63 |
| black glass ⊛ GGX 0.019 | 1.52 | 1.59 | 1.67 | 1.85 | 2.11 | 2.71 | 3.35 |
| black glass, slopes scaled ×1.38 | 1.51 | 1.59 | 1.78 | 2.06 | 2.49 | 3.40 | 4.29 |
| Aluminum, measured | 1.74 | 1.80 | 1.92 | 2.10 | 2.48 | 3.45 | 4.22 |
| black glass ⊛ GGX 0.027 | 1.73 | 1.80 | 1.86 | 2.03 | 2.29 | 2.85 | 3.48 |
| black glass, slopes scaled ×1.56 | 1.70 | 1.80 | 2.01 | 2.32 | 2.81 | 3.84 | 4.84 |

Scaling the slopes instead (P(m/s)/s², the other natural way to widen a lobe) also scales the instrument's response, and
overshoots from 40° on. The convolution also keeps the in-plane/cross-plane difference (aluminium at θ_d 50°: measured
2.43° / 1.97°, glass ⊛ kernel 2.19° / 2.01°, scaled glass 2.77° / 2.31°). Aluminium widens more than the glass at 70°–75°,
which a kernel independent of θ_d does not reproduce. The tails do not follow GGX: aluminium has 6.6 × 10⁻³ of its peak at
θ_h = 10° and 1.4 × 10⁻³ at 20°, the glass ⊛ GGX 2.0 × 10⁻³ and 2.1 × 10⁻⁴. This haze belongs to the material, and
`gloss` keeps it because it convolves the material itself.

Consequences:

- The roughness of the surfaces themselves is smaller than alpha_eq suggests (aluminium: a GGX kernel of 0.027 on top of the
  instrument). From about alpha_eq 0.1 (hwhm 3.7°) on, the instrument's part of the width is small.
- A convolution keeps the measurement's look, and widens a lobe consistently at every θ_d: at alpha_eq 0.21 the hwhm is
  7.9°–8.5° from θ_d 20° to 70°, nearly constant as for a surface.
- Nothing narrower than the glass lobe was measured. `gloss` refuses a target below the material's width; a smoother
  version of a material is a swap of its optics onto a sharper measured shape (for example aluminium's colour,
  reflectance and polarization on the lobe of Fake gold, `brdfkit swap`).

## The knob against the MERL components

`brdfkit edit --pc 1:T` moved the lobe along the first MERL component ([LATENT_EDIT.md](LATENT_EDIT.md)). On the same clean
aluminium (PC edits with the albedo kept at 45°, as in the demo; `gloss` without rescaling):

| Lobe | hwhm at θ_d 20° / 30° / 50° / 70° | D at θ_h = 0 / its maximum | D / peak at θ_h 10° / 20° / 40° | Albedo at θ_i 0° / 45° / 75° |
|---|---|---|---|---|
| aluminium | 1.74 / 1.80 / 2.10 / 3.45 | 1 | 6.6e-3 / 1.4e-3 / 2.0e-4 | 0.77 / 0.90 / 0.70 |
| PC1 −1 | 0.21 / 0.23 / 0.26 / 0.37 | 1 | 2.0e-4 / 2.9e-5 / 2.8e-6 | 0.64 / 0.90 / 0.63 |
| PC1 +0.5 | 2.37 / 2.47 / 2.84 / 4.51 | 0.53 | 2.0e-2 / 5.3e-3 / 8.8e-4 | 0.86 / 0.90 / 0.77 |
| PC1 +1 | 2.84 / 2.98 / 3.47 / 5.09 | 0.21 | 4.5e-2 / 1.5e-2 / 2.9e-3 | 0.96 / 0.90 / 0.87 |
| gloss, alpha_eq 0.08 | 2.92 / 2.96 / 3.18 / 4.15 | 1 | 3.3e-2 / 5.4e-3 / 6.9e-4 | 0.76 / 0.88 / 0.67 |
| gloss, alpha_eq 0.13 | 4.82 / 4.84 / 5.01 / 5.68 | 1 | 1.4e-1 / 2.3e-2 / 2.8e-3 | 0.76 / 0.86 / 0.65 |
| gloss, alpha_eq 0.21 | 7.95 / 7.95 / 8.04 / 8.47 | 1 | 3.7e-1 / 9.2e-2 / 1.3e-2 | 0.76 / 0.80 / 0.64 |
| gloss, alpha_eq 0.34 | 13.5 / 13.5 / 13.4 / 13.6 | 1 | 6.6e-1 / 2.9e-1 / 6.3e-2 | 0.72 / 0.70 / 0.62 |
| gloss, alpha_eq 0.55 | 25.2 / 25.1 / 24.5 / 21.0 | 1 | 8.7e-1 / 6.3e-1 / 2.6e-1 | 0.60 / 0.55 / 0.57 |

(hwhm in degrees, from the maximum of D outward, for the gloss lobes with the applied shadowing divided out; at θ_d 70°
the widest lobe is cut by the horizon. Fake gold along the knob has the same widths within 0.05° up to θ_d 50° and
alpha_eq 0.34.)

What PC1 does to these lobes:

- **Widening makes a ring.** Along +PC1 the centre of the lobe drops below its shoulders: at PC1 +1, D at θ_h = 0 is 0.21 of
  its maximum at θ_h ≈ 1.6°. Under a point light the highlight is a ring (`out/gloss/gloss_axis.png`, panel e). The change
  that the MERL basis applies near the peak is meant for MERL's resolved peaks; on the flat, wider hpBRDF peak it digs a hole.
  `gloss` keeps the maximum at the centre.
- **It changes more than the width.** At about the same hwhm (PC1 +1 and alpha_eq 0.08, both 3.0°), PC1 has 3–4 times
  more light at θ_h 20°–40°, and its albedo at normal incidence rises from 0.77 to 0.96 while the albedo at 45° is held:
  it adds a broad, paint-like component, which is why the nearest MERL materials become paints.
- **It keeps the instrument's θ_d trend** (×1.7 from θ_d 30° to 70° at PC1 +1), where a rougher surface has none.
- **It reaches little.** PC1 +1 is 3.0° wide; larger shifts leave the plausible materials. Toward −PC1 it narrows the lobe
  to 0.2°, five times below the instrument's resolution: a sharpened flat top, not a measured shape.
- PC3, read from MERL as the lobe width, narrows the aluminium toward +PC3 (0.99° at +0.5, 0.38° at +1) and widens it
  toward −PC3, also into a ring (−1.5: 2.55°, centre 0.26 of the maximum).

The knob reaches a matte metal (alpha_eq 0.55, hwhm 25°) with one parameter that means the same for every material.

## Rendered

GoLD_rayt (spectral build) under an environment light (`out/gloss/row_*.json`, contact sheet `out/gloss/gloss_sheet.png`),
the renderer's importance sampling set to max(0.1, alpha_eq) (`alphaSample`, which changes only the noise). Per
sphere: the mean luminance, the reflection contrast (95th over 5th percentile of the luminance on the sphere), the
brightest 0.5% of the luminance, and the degree of linear polarization (DoLP, luminance weighted).

Aluminium along the knob:

| alpha_eq (hwhm) | Luminance | Contrast | Brightest 0.5% | DoLP |
|---|---|---|---|---|
| aluminium, 0.049 (1.8°) | 0.253 | 23 | 6.6 | 0.034 |
| 0.08 (3.0°) | 0.245 | 25 | 5.1 | 0.036 |
| 0.13 (4.8°) | 0.243 | 21 | 3.6 | 0.034 |
| 0.21 (8.0°) | 0.230 | 15 | 2.1 | 0.032 |
| 0.34 (13.5°) | 0.202 | 11 | 1.1 | 0.029 |
| 0.55 (25°) | 0.150 | 8.2 | 0.54 | 0.025 |

Fake gold along the knob:

| alpha_eq (hwhm) | Luminance | Contrast | Brightest 0.5% | DoLP |
|---|---|---|---|---|
| Fake gold, 0.031 (1.1°) | 0.176 | 29 | 5.1 | 0.058 |
| 0.08 (3.0°) | 0.171 | 36 | 4.1 | 0.056 |
| 0.13 (4.8°) | 0.171 | 27 | 2.8 | 0.053 |
| 0.21 (8.0°) | 0.163 | 17 | 1.5 | 0.046 |
| 0.34 (13.5°) | 0.144 | 13 | 0.78 | 0.039 |
| 0.55 (25°) | 0.108 | 8.7 | 0.36 | 0.035 |

Both spheres go from a polished to a matte metal; the colour stays (aluminium RGB / luminance 1.15–1.18, 0.98, 0.72;
Fake gold 1.57–1.60, 0.90–0.91, 0.22–0.23). The same alpha_eq gives the same highlight size on both materials: the
brightest 0.5% falls to about 0.3 of its smooth value at 0.21 on each (aluminium 6.6 → 2.1 with its own 1.8° start, Fake gold
5.1 → 1.5 from 1.1°). The highlights spread first while the dark reflections stay dark (the contrast even rises a little at 3°,
where the blurred highlights still cover little of the sphere); from about 8° on the reflected room blurs away and the
contrast falls. The luminance falls by 40% at alpha_eq 0.55, the single-scattering loss above. The DoLP falls a little although N is unchanged: a wider lobe mixes, in each pixel, light
reflected with different planes of incidence.

MERL PC1 against the knob at the same width:

| Lobe (hwhm) | Luminance | Contrast | Brightest 0.5% | DoLP |
|---|---|---|---|---|
| PC1 +0.5 (2.5°) | 0.256 | 18 | 5.2 | 0.036 |
| gloss (2.5°) | 0.248 | 23 | 5.9 | 0.042 |
| PC1 +1 (3.0°) | 0.259 | 13 | 4.4 | 0.039 |
| gloss, alpha_eq 0.08 (3.0°) | 0.251 | 24 | 5.3 | 0.048 |

At the same width the PC1 spheres look milky: their contrast is that of the knob at about 6° (PC1 +0.5) and 9° (+1),
while their highlights are as small as at 2.5°–3°. This is the broad component PC1 adds (table above): a veil over the
whole sphere, not a blur of the reflections.

## Limits

- **Polarization does not follow the roughness.** N is kept. A rough surface depolarizes more, mostly at small θ_d, through
  multiple scattering (hpBRDF paper, Fig. 12); a roughened material keeps the polarization of its smooth original.
- **Single scattering.** A rough metal returns part of the light that single scattering loses (multiple scattering between
  microfacets); `gloss` does not, so the albedo of the widest lobes is lower than a real surface's: at alpha_eq 0.21 the albedo
  at 45° drops by 11% (0.90 → 0.80), by 39% at 0.55 (Fake gold: 13% and 40%). `--keep-albedo` holds it at one angle.
- **The instrument's response stays.** It is about 1.1° at θ_d 20°–30° and widens with θ_d; for rough lobes it is negligible,
  for nearly smooth ones it is most of the width.
- **Widening only.** A lobe narrower than the material's cannot be recovered from the measurement (see above).
- **Perceptual scale.** alpha_eq is a physical scale. Equal steps in log alpha_eq are a reasonable start for stimuli, but
  whether they look equally spaced has to be measured (for example with maximum-likelihood difference scaling).
