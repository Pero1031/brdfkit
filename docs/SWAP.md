# Material swap: the shape of one material, the colour, polarization and reflectance of another

Implementation in `brdfkit/analysis/swap.py`.

```bash
# clean models first: calibrated on black glass and symmetrized over phi_d (edit without a shift)
python -m brdfkit edit out/factor/Aluminum_s0.npz  --model out/demo/merl_pca/merl_pca.npz --out out/swap/Aluminum_clean.npz  --calib out/calib/black_glass_cal.npz
python -m brdfkit edit out/factor/Fake_gold_s0.npz --model out/demo/merl_pca/merl_pca.npz --out out/swap/Fake_gold_clean.npz --calib out/calib/black_glass_cal.npz
# aluminium's lobe with everything optical of Fake gold, or only its colour
python -m brdfkit swap out/swap/Aluminum_clean.npz --source out/swap/Fake_gold_clean.npz --out out/swap/al_as_Fake_gold.npz
python -m brdfkit swap out/swap/Aluminum_clean.npz --source out/swap/Fake_gold_clean.npz --parts S --out out/swap/al_gold_S.npz
# or Fresnel of an n, k table as the source (needs a calibrated shape)
python -m brdfkit swap out/swap/Aluminum_clean.npz --nk data/kaist/refidx/Al_Rakic.yml --film-nm 20 --out out/swap/al_physical.npz
```

The result is an ordinary factor model: `compress` and `export` bake it for GoLD_rayt, and `edit` can move its lobe on
top of the swap (for example a rougher gold).

## Why this split

A specular (microfacet) reflection is f = F(θ_d, λ) · D G / (4 cos θ_i cos θ_o). The Fresnel factor F depends only on
the angle of incidence on the microfacet, θ_d, and on the optical constants: it gives the reflectance, its colour and its
polarization. The rest depends on the surface geometry (roughness). The factor model ([FACTOR.md](FACTOR.md)) has the
same structure: m00 = A × S with S normalized to luminance 1, and m_ij / m00 = N(θ_d, λ). So A = level(θ_d) × shape,
where level is the luminance reflectance along θ_d. A swap keeps the shape of material X and takes any of these from
material Y:

| Part | What it is | Swap |
|---|---|---|
| `S` | colour (68 bands, luminance 1) | S' = S_Y: the luminance stays |
| `N` | polarization (15 normalized Mueller elements) | N' = N_Y |
| `level` | reflectance along θ_d | A' = A_X · level_Y / level_X |

The level is the lobe energy around the specular direction (θ_h < 30°, smoothed over 3° of θ_d), the quantity that the
calibration compares with Fresnel ([CALIBRATION.md](CALIBRATION.md)). The instrument error K(θ_d) is the same for
every hpBRDF table, so it cancels in level_Y / level_X; the two models only have to be in the same calibration state,
and `swap` refuses a mix. A Fresnel source (`--nk`) gives the true reflectance, so it needs a calibrated shape.

Under unpolarized light and one bounce to the camera, N does not change the image's radiance, only its polarization.
It matters when a polarizer is in the view, and after a second reflection (a sphere reflected in glass or water).

## The four downloaded materials

From their own measurements (factor models, calibrated):

| Material | Type | Retardance at 550 nm | Reflectance at θ_d 30° / 50° / 70° | Colour: S at 450 / 550 / 650 nm |
|---|---|---|---|---|
| Aluminum | metal | smooth fall, principal angle 73.5° | 1.03 / 0.89 / 0.82 | 0.86 / 1.04 / 0.73 |
| Fake gold | metal | smooth fall, principal angle 70.6° | 0.58 / 0.61 / 0.68 | 0.34 / 0.96 / 1.30 (yellow) |
| SUJ2 (bearing steel) | metal | smooth fall, principal angle 66° | 0.46 / 0.51 / 0.72 | 0.94 / 0.90 / 0.94 |
| Black glass | dielectric | 146° step at Brewster (57.5°), D = 1 there | 0.045 / 0.062 / 0.19 | about flat |

## Result: the optics move, the shape stays

GoLD_rayt (spectral build) under an environment light, all models calibrated and symmetrized (`out/swap/row_*.json`).
For each sphere: the mean luminance, the mean RGB divided by the luminance, and the luminance-weighted degree of linear
polarization (DoLP):

| Sphere | Luminance | RGB / luminance | DoLP |
|---|---|---|---|
| Fake gold, measured | 0.164 | 1.61 0.90 0.24 | 0.068 |
| aluminium shape + Fake gold | 0.160 | 1.60 0.90 0.25 | 0.059 |
| SUJ2, measured | 0.108 | 1.41 0.91 0.68 | 0.127 |
| aluminium shape + SUJ2 | 0.109 | 1.41 0.91 0.68 | 0.124 |
| black glass, measured | 0.0136 | 1.47 0.88 0.76 | 0.586 |
| aluminium shape + black glass | 0.0134 | 1.47 0.88 0.75 | 0.581 |

The swapped spheres match the measured ones within 3% in luminance and colour, and within 15% in DoLP. What still
differs is visible: Fake gold and black glass are smoother than aluminium (lobe half width 1.3° against 1.9°), so their
measured spheres show sharper reflections; SUJ2 has the width of aluminium, and its swapped and measured spheres look
alike. The shape is the one thing the swap keeps on purpose.

With a Fresnel source (`--nk Al_Rakic.yml --film-nm 20`) on the calibrated aluminium, the level changes by a factor of
1.06: the calibrated measurement is within 6% of Fresnel for aluminium with a 20 nm oxide film.

One part at a time (aluminium shape, parts of Fake gold):

| Sphere | Luminance | RGB / luminance | DoLP |
|---|---|---|---|
| aluminium | 0.235 | 1.16 0.98 0.70 | 0.039 |
| + colour (S) | 0.235 | 1.61 0.90 0.24 | 0.043 |
| + colour, reflectance (S, level) | 0.161 | 1.59 0.90 0.24 | 0.049 |
| + colour, reflectance, polarization (S, level, N) | 0.162 | 1.58 0.91 0.24 | 0.064 |

- The colour swap changes the hue and keeps the luminance, by construction. The sphere looks gold, but brighter than the
  measured gold.
- The reflectance swap brings the luminance down to that of Fake gold (0.235 → 0.161).
- The polarization swap leaves the image unchanged (0.161 → 0.162, within the render noise). It changes the DoLP from
  0.049 to 0.064, and the DoLP image changes from aluminium's bright ring at the rim to the pattern of Fake gold
  (`out/swap/swap_sheet.png`, last row).

The three parts can therefore be set independently: colour and brightness change what an observer sees, and polarization
only what a polarization camera, a polarizing filter or a second reflection sees.

## Limits

- The level is the lobe energy, which holds the reflection of metals and glass. Materials with a strong diffuse part
  (the rough plastics of hpBRDF) need a specular and a diffuse component; a swap would mix them.
- S and N carry the measurement errors of their source ([CALIBRATION.md](CALIBRATION.md): about 10%–25% in the spectrum
  and up to about 0.18 in diattenuation at low θ_d); a Fresnel source (`--nk`) avoids them, at the cost of the measured
  detail (for example the oxide film of the real aluminium).
- A roughened or sharpened material keeps its source's polarization. In reality roughness also changes the
  depolarization ([LATENT_EDIT.md](LATENT_EDIT.md), Limits).
