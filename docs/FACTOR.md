# Factorized hpBRDF: grey angular shape × spectrum × polarization

For material editing in a spectral, polarized renderer it helps to split a measured hpBRDF into parts that can be edited
separately: the angular shape of the lobe (which a MERL-like basis covers well), the spectrum, and the polarization.
`brdfkit factor` fits

    m00(φ_d, θ_d, θ_h, λ) = A(φ_d, θ_d, θ_h) · S(θ_d, λ)
    m_ij / m00            = N_ij(θ_d, λ)        (the 15 other Mueller elements)

Implementation in `brdfkit/analysis/factor.py`.

```bash
python -m brdfkit factor   data/kaist/hpbrdf/Aluminum.hpbrdf --out out/factor/Aluminum.npz --smooth 0,5 --floor out/floor/Black_glass_floor.npz
python -m brdfkit compress out/factor/Aluminum_s0.npz out/factor/Aluminum_factor_pca.npz      # a model loads like a table
python -m brdfkit export   out/factor/Aluminum_factor_pca.npz out/factor/Aluminum_factor.hpbpca
```

## Why this form

The table's Mueller matrices use Stokes frames with the x reference vectors in the plane of reflection (Mitsuba
`measured_polarized` convention). For a mirror microfacet that plane is its plane of incidence and the angle of incidence on
it is θ_d, so a single-bounce specular reflection is a Fresnel Mueller matrix in (θ_d, λ) times a grey distribution term in
the angles. A is in luminance units: S is normalized to CIE Y = 1 for every θ_d, so A can be compared with, or replaced by,
a MERL lobe.

## Fit

Two streaming passes over the table (floor removed first, see [FLOOR.md](FLOOR.md)):

1. Per θ_d, the Gram matrix of the m00 spectra over all bins; S is its top eigenvector, the energy-optimal rank-1 spectrum.
   N is the m00-weighted mean of clip(m_ij / m00, −1, 1) over φ_d and θ_h.
   Outside the measured θ_d range (20°–80°; the rest of the hpBRDF table is inpainted, see [LU_CHIPMAN.md](LU_CHIPMAN.md))
   S and N are copied from the nearest measured θ_d. `--smooth` optionally smooths S along θ_d (Gaussian, degrees).
2. A = m00 · S / (S · S) per bin, and the reconstruction errors.

Errors are energy weighted: for m00 the relative RMS sqrt(Σ (m00 − A S)² / Σ m00²), for N the RMS over the 15 elements
weighted by m00. The model is 6.8 MB (A in float16), against 13 GB for the table and 231 MB for the PCA compression.

## Results (floor removed)

m00 relative RMS error:

| | θ_h < 5° | 5°–20° | ≥ 20° | all, measured θ_d | all θ_d |
|---|---|---|---|---|---|
| Aluminum, S per θ_d | 0.060 | 0.26 | 0.32 | **0.062** | 0.114 |
| Aluminum, S smoothed 5° | 0.096 | 0.26 | 0.31 | 0.097 | 0.132 |
| Black glass, S per θ_d | 0.063 | 0.47 | 0.60 | **0.077** | 0.130 |
| Black glass, S smoothed 5° | 0.077 | 0.47 | 0.59 | 0.089 | 0.137 |

N RMS error (same for both S variants):

| | θ_h < 5° | 5°–20° | ≥ 20° | all, measured θ_d | all θ_d |
|---|---|---|---|---|---|
| Aluminum | 0.010 | 0.056 | 0.21 | 0.016 | 0.030 |
| Black glass | 0.008 | 0.039 | 0.40 | 0.015 | 0.041 |

- Near the specular ridge, which carries almost all the energy, the model is within 6% in m00 and 0.01 in the normalized
  Mueller elements (whose magnitude is about 0.4). For comparison the 8 + 32 component PCA compression has about 10% m00
  error and 0.011–0.012 polarization error near the ridge ([PCA.md](PCA.md), [FLOOR.md](FLOOR.md)).
- The shoulders and the far field (θ_h ≥ 5°) do not follow the model: their spectrum and polarization differ from the
  ridge (multiple scattering, depolarization). They are dark, so they barely count in the energy-weighted total.
- The "all θ_d" column includes the inpainted θ_d ranges, where the model copies the edge values while the table holds the
  inpainting; neither is a measurement there.
- Smoothing S along θ_d costs 1–3.5 points of m00 error, so the irregular θ_d dependence of S is consistent in the data
  rather than noise of the fit. It does not look physical (on black glass the 798 nm / 550 nm ratio of S goes 1.17, 0.94,
  1.15, 1.53, 1.56 from θ_d 20° to 79°, where glass is nearly flat); a calibration difference between illumination angles
  is the likely cause, not yet confirmed.

## Use for editing

A, S and N are independent arrays, so a material can be recombined: the S and N of one material (or of a physical model,
n, k and a film) with the A of another, including a MERL lobe on the hpBRDF grid. The θ_h ≥ 5° part is the weakest; a
second, depolarized component would be the next refinement if it shows in renderings.
