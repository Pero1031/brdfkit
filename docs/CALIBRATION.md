# Black glass as a built-in reference: the radiometric error of the hpBRDF tables

A polished dielectric reflects exactly as Fresnel predicts from one real index n(λ), and the hpBRDF set contains one:
black glass. Comparing it with Fresnel shows how far the tables are from a correctly scaled, correctly polarized BRDF.
Comparing a second material with a known index (aluminium) separates errors of the instrument (common to both) from
properties of the samples. Implementation in `brdfkit/analysis/calibrate.py`.

```bash
python -m brdfkit factor    data/kaist/hpbrdf/Black_glass.hpbrdf --out out/factor/Black_glass.npz --floor out/floor/Black_glass_floor.npz
python -m brdfkit calibrate out/factor/Black_glass_s0.npz --out out/calib/black_glass_cal.npz \
                            --check out/factor/Aluminum_s0.npz --refidx data/kaist/refidx/Al_Rakic.yml --film-nm 20
python -m brdfkit edit      out/factor/Aluminum_s0.npz --model out/merl_analysis/log/merl_pca.npz --out out/calib/Aluminum_sym_cal.npz \
                            --calib out/calib/black_glass_cal.npz          # no shift: the cleaned, calibrated material
```

`calibrate` writes the calibration (.npz), a report (.txt) and a figure (.png). `edit --calib` divides the lobe A by
K(θ_d) before the symmetrization and the edit. `--calib-keep-level 30` corrects only the θ_d trend and keeps the level at
θ_d = 30°. `examples/latent_demo.py --calib` does the same for the demo.

## Method

All three checks use the factor model of the glass (m00 = A × S, m_ij/m00 = N, [FACTOR.md](FACTOR.md)), with the floor
removed. They use only the measured θ_d range (20°–79°).

1. **Polarization.** N(θ_d, λ) is rotated into the (s, p) bases at the specular geometry. It is fitted with the Fresnel
   Mueller matrix of a dielectric (m01, m10, m22, m33), which gives n(λ). Lu–Chipman D and R of the data are then compared
   with those of the fit.
2. **Reflectance.** The energy of the lobe around the specular direction is
   Q(θ_d) = 2π · 4 cos²θ_d · ∫ Ā(θ_d, θ_h) cos θ_h sin θ_h dθ_h, taken over θ_h < 10° with Ā the mean of A over φ_d.
   For a narrow lobe of a correctly scaled BRDF, this is the directional albedo at incidence θ_d. For the glass it should
   equal the Fresnel reflectance F_Y(θ_d) of the fitted index. K = Q / F_Y is the radiometric error of the table.
   In the hpBRDF instrument the camera is fixed and θ_d is half the light–camera angle, so one θ_d is one illumination
   angle. K is smoothed over 3° of θ_d and copied from the edges outside the measured range.
3. **Spectral shape.** S(θ_d, λ) is compared with F(θ_d, λ) / F_Y. Glass reflects almost the same at every wavelength.

## Result 1: the lobe energy falls far below Fresnel, the same way for glass and aluminium

The directional albedo divided by the Fresnel reflectance should be 1. It was computed by Monte Carlo over the hemisphere
(`latent_edit.albedo`, an exact integral, not the approximation above), for black glass (fitted n) and for aluminium
(Rakić n, k with a 20 nm oxide film, [LU_CHIPMAN.md](LU_CHIPMAN.md)):

| Incidence | 25° | 35° | 45° | 55° | 65° | 75° |
|---|---|---|---|---|---|---|
| Black glass: albedo / Fresnel | 0.46 | 0.37 | 0.27 | 0.19 | 0.14 | 0.09 |
| Aluminium: albedo / Fresnel | 0.43 | 0.37 | 0.27 | 0.19 | 0.13 | 0.10 |
| Aluminium / glass | 0.92 | 0.98 | 0.99 | 0.98 | 0.91 | 1.10 |

The two materials differ by a factor of 20 in reflectance (4% against 90%), in lobe width (half width 1.3° against 1.9°)
and in the shape of their shoulders. Their errors agree within 10%, so the error belongs to the tables, not to the samples. The tables carry about half
the Fresnel energy at 20°–25° and a tenth at 75°. Up to 45° the error follows 0.55 cos²θ_d closely, and it is flatter
beyond. K estimated from the glass:

| θ_d | 20° | 25° | 30° | 35° | 40° | 45° | 50° | 55° | 60° | 65° | 70° | 75° |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K | 0.47 | 0.43 | 0.38 | 0.35 | 0.30 | 0.25 | 0.21 | 0.17 | 0.15 | 0.12 | 0.10 | 0.08 |

Correcting the aluminium with the glass's K brings its albedo / Fresnel from 0.43–0.10 to **1.01, 1.05, 1.00, 0.92,
0.85, 0.86** at 25°–75°. This check is independent, because the aluminium was not used to estimate K. At 65°–75° the
remaining shortfall is about 15%. Energy that would leave below the horizon is lost there, θ_d > 79° is inpainted, and
the integral over θ_h < 10° depends on the cut-off near grazing (5° instead of 10° lowers K at 75° by 25%).

The cause is not known. Candidates are an illumination-angle-dependent radiometric calibration, or tables that already
contain cosine factors. The Mitsuba 3 `measured_polarized` plugin, whose file format the hpBRDF release uses, multiplies
the table value by cos θ_o, so it treats the table as a plain BRDF. GoLD_rayt and brdfkit do the same. Rendered as
delivered, hpBRDF materials are therefore about twice too dark and lose the brightening toward grazing that every smooth
surface has. Metal and glass spheres get dark rims.

The renders below are from GoLD_rayt (spectral build), with spheres under an environment light, both versions with the
lobe symmetrized (`out/calib/compare_env.json`). Rim / centre is the mean luminance at 0.85–0.97 of the radius over
that within 0.5, upper half only. That part reflects the room, not the floor:

| | Mean luminance | Rim / centre |
|---|---|---|
| Aluminium, as delivered | 0.073 | 0.26 |
| Aluminium, calibrated | 0.213 | 0.45 |
| Black glass, as delivered | 0.0037 | 0.80 |
| Black glass, calibrated | 0.0131 | 1.77 |

Calibrated, the aluminium looks like polished metal rather than dark gunmetal. The black glass reflects the room and
brightens toward its rim, as a glass ball does.

## Result 2: the spectral errors of the glass do not repeat on aluminium

S / (Fresnel shape) of the glass, at a few wavelengths (1 = correct):

| θ_d | 414 nm | 462 nm | 550 nm | 702 nm | 798 nm | 942 nm |
|---|---|---|---|---|---|---|
| 20° | 1.01 | 1.02 | 0.76 | 0.84 | 0.84 | 0.61 |
| 30° | 1.18 | 0.88 | 1.02 | 0.75 | 0.83 | 0.70 |
| 40° | 1.07 | 1.01 | 0.93 | 0.84 | 0.86 | 0.68 |
| 50° | 1.14 | 1.02 | 0.91 | 0.92 | 1.00 | 0.73 |
| 60° | 0.90 | 1.04 | 0.93 | 0.93 | 1.15 | 0.81 |
| 70° | 0.99 | 1.08 | 0.89 | 0.90 | 1.19 | 0.88 |

The glass's spectrum is off by 10%–25% in the visible and by up to 40% at the near-infrared edge, irregularly along θ_d.
The aluminium's deviations from its own Fresnel shape have the same size but do not follow the glass's: at 30° and
702 nm the glass reads 0.75 and the aluminium 1.25. A spectral correction taken from the glass would therefore not fix
other materials, and `apply` leaves S alone. Smoothing S along θ_d ([FACTOR.md](FACTOR.md)) or replacing it with a
physical model is the better route.

## Result 3: polarimetric accuracy at the specular

- **Index.** The glass's polarization gives n = 1.503–1.555 over wavelength (median 1.529), plausible for glass. It
  scatters by ±0.015 from band to band without the expected dispersion. The fit residual (RMS of m01, m10, m22, m33) is
  0.04–0.11.
- **Diattenuation.** At 550 nm, D is off by −0.18 to +0.05 at θ_d 20°–30°: the glass reads about 0.2 where Fresnel rises
  from 0.17 to 0.39. The error is at most 0.09 at 60°–79°. Over all wavelengths, the median error is 0.055 at 20°–30° and
  0.018 at 60°–79°.
- **Retardance.** The error is at most 2° at 20°–50° and up to 9° at 65°–79°.
- **m02**, which an ideal reflector has at 0, drifts from +0.1 at 25° to −0.3 at 55°–70°. The aluminium does not show
  the same drift (|m02| ≤ 0.08), so it is not a single rotation of the polarizer frames.

For the aluminium analysis ([LU_CHIPMAN.md](LU_CHIPMAN.md)), this means the 20 nm oxide film, read from the retardance,
holds. The excess diattenuation at 20°–30° is within the instrument's error. The excess at 60°–79° is beyond it.

## Limits

- K comes from one sample. Black glass is the only dielectric of known physics in the set, and its absolute level rests on
  its fitted index (n ≈ 1.53, F_Y = 0.043 at 20°).
- The correction depends on θ_d only. That is the right form for an error per illumination angle. If the tables carry
  cosine factors of θ_i and θ_o instead, the shoulders of a broad lobe (θ_h ≫ 1°) are corrected only approximately.
- Beyond 65° the estimate depends on the θ_h cut-off, and beyond 79° it is copied from the edge.
