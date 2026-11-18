# Polarimetric analysis of measured aluminium (hpBRDF) by Lu–Chipman decomposition

Results of `python -m brdfkit lu-chipman out/pca/Aluminum_pca.npz --refidx data/kaist/refidx/Al_Rakic.yml data/kaist/refidx/Al_McPeak.yml --out out/analysis --raw data/kaist/hpbrdf/Aluminum.hpbrdf`.
Implementation in `brdfkit/analysis/lu_chipman.py`; output in `out/analysis/` (3 figures, `lu_chipman_specular.npz`, `lu_chipman_report.txt`).

## Method

The Mueller matrix is factored as M = M_Δ · M_R · M_D (depolarizer × retarder × diattenuator; Lu & Chipman 1996), giving

- diattenuation D (0–1): the difference between the s and p reflectances, (|r_s|² − |r_p|²)/(|r_s|² + |r_p|²)
- retardance R (0–180°): the phase shift between s and p. An ideal reflection has 180° at normal incidence, 90° at the principal angle, 0° at grazing
- depolarization Δ (0–1): the fraction of polarization lost. 0 for an ideal Fresnel surface
- and the Gil–Bernabeu depolarization index DI (1 = polarization preserved)

applied to the Mueller matrices decoded from the PCA-compressed table (rotated into the same (s, p) bases as GoLD_rayt).
The references are the Fresnel values from the complex refractive index of Al (Rakić 1995, McPeak 2015) and a thin-film
model with one dielectric layer of index 1.65 on Al (the native oxide Al₂O₃); the film thickness is fitted to D(θ) and R(θ).

**Measured range.** The hpBRDF light source is rotated through 40°–160° in 5° steps (Moon et al. 2025, Sec. 4), so only
θ_d ≈ 20°–80° is measured; the rest of the table, notably low θ_d, is filled by 3D Gaussian-convolution inpainting.
In the specular configuration used below θ = θ_d, and the decoded D, R, Δ are constant to within 0.003 for θ ≤ 12.5°:
those values are copies of the measured edge, not measurements. This document therefore draws no conclusion from θ < 20°.
(The film-thickness fit below uses θ = 5°–75°, so the inpainted points enter it as repeats of the edge value; the
principal angles that carry the result lie at 72°–77°, inside the measured range.)

Validation: the theoretical Fresnel matrices give Δ = 0 and D, R equal to the analytic values, and 50 matrices synthesized
from known parameters are recovered. Decomposing the original 13 GB table directly (the θ_h = 0 slice) and decomposing the
decoded values agree within 0.02 in D, 1° in R and 0.01 in Δ.

## Result 1: the retardance R is explained by "Al + about 20 nm of dielectric film"

Principal angle (where R crosses 90°) in the specular configuration (θ_i = θ_o):

| Wavelength | Measured | Bare Al (Rakić) | Bare Al (McPeak) | Al + 20.2 nm film |
|---|---|---|---|---|
| 446 nm | 72.2° | 79.9° | 78.5° | 71.1° |
| 550 nm | 73.6° | 81.6° | 80.6° | 74.2° |
| 646 nm | 72.5° | 82.8° | 81.9° | 76.2° |
| 798 nm | 77.1° | 83.5° | 82.8° | 78.2° |

Bare Al puts the principal angle 8–10° higher than measured; one film layer brings it within 1–4°.
The whole R(θ) curve follows: the dashed thin-film model nearly overlaps the solid measured curve (`lu_chipman_specular.png`, middle).
Fitting the thickness per wavelength gives a monotonic increase from 14.5 nm (414 nm) to 29 nm (950 nm), so the assumption
"one layer with a constant index of 1.65" is crude and the film's dispersion or roughness also contributes; still,
**as far as the retardance is concerned, the measured aluminium behaves like Al with an oxide film**.
The black glass ([CALIBRATION.md](CALIBRATION.md)) shows a retardance error of at most 2° at 20°–50° and up to 9° at
65°–79°. That moves the aluminium's principal angle by about 1–2°, far less than the measured 8–10° shift, so this
conclusion holds.
R is smooth along both wavelength and angle and is the most reliable quantity in this analysis.

## Result 2: the diattenuation D is explained by neither Fresnel nor the film

| θ | measured D, 550 nm | bare Al | 20 nm film | measured D, 800 nm | bare Al |
|---|---|---|---|---|---|
| 20° | 0.091 | 0.006 | 0.006 | 0.239 | 0.009 |
| 25° | 0.104 | 0.009 | 0.009 | 0.199 | 0.014 |
| 45° | 0.187 | 0.031 | 0.031 | 0.252 | 0.050 |
| 70° | 0.321 | 0.099 | 0.073 | 0.432 | 0.167 |

The measured D is 3 to almost 30 times the theory. At the smallest measured angles (θ = 20°–30°) it is 0.10 (median over all
wavelengths; range 0.02–0.47), where Fresnel and the film give at most 0.02.
For an isotropic surface D → 0 towards normal incidence, and at 20° it is already below 0.01 for Al with or without the film;
the excess is therefore not explained by the material's optical constants but indicates either

- an anisotropic sample (rolled or brushed aluminium reflects differently along and across the grooves, giving diattenuation even near normal incidence), or
- a systematic error of the polarimeter (averaged out by the dataset's isotropy assumption).

Earlier versions of this document quoted D at 10°; that angle lies in the inpainted range (see Method) and only repeats
the value at the measured edge. Whether D goes to 0 at normal incidence cannot be tested with this dataset.

**Checked against black glass** ([CALIBRATION.md](CALIBRATION.md)). The hpBRDF black glass should follow Fresnel of a
dielectric exactly. Its D at 550 nm is off by −0.18 to +0.05 at θ = 20°–30°, and by at most 0.09 at 60°–79°. So the
aluminium's excess of about 0.09 at 20°–30° is within the instrument's error there and is not established as a property
of the sample. Its excess of 0.22–0.30 at 60°–79° is larger than the glass error. Later checks on six materials
([QUIRKS.md](QUIRKS.md)) withdraw it too: three dielectrics follow Fresnel, but the metals' diattenuation does not, in
both directions (aluminium five times too high, SUJ2 eight times too low against iron), and differs from their
polarizance. The excess is a problem of the metal measurements, not a property of this sample.

D is also noisy along wavelength (`lu_chipman_spectral.png`, left).
The film barely changes D (the dashed film curve is almost identical to the dotted bare one), so **the thickness information comes from R alone**.

## Result 3: the depolarization Δ is small and grows away from the specular ridge

- In the specular configuration Δ = 0.00–0.17, DI = 0.83–1.00: a reasonable range for a metal of low roughness.
- Δ grows with angle (550 nm: 0.05 at 20°, 0.17 at 70°).
- Off the specular ridge Δ rises to 0.25–0.5 and exceeds 0.7 where the light is near the normal and the viewer near grazing
  (`lu_chipman_map550.png`, right): multiple scattering and the tail component destroy polarization.
- Around 800 nm Δ goes negative (−0.04) and DI exceeds 1. The matrix slightly leaves the physically admissible range
  (sum of squares of |m_ij| > 3 m00²), presumably from measurement noise plus compression error.
  This argues for a physically constrained compression (projection onto passive matrices).

## Note: symmetry of the table

At θ_h = 0, isotropy means the Mueller matrix should not depend on φ_d, yet the φ_d = 0 and φ_d = π slices differ by up to
11% at 20° (the edge of the measured range) and 1% at 60°. Off the specular ridge, R is also asymmetric in θ_i and θ_o. Both show that the table is not
strictly isotropic or reciprocal; harmless for rendering, but analyses should average over them.

## Summary

| Quantity | Measured vs ideal Al | Interpretation |
|---|---|---|
| Retardance R | principal angle 8–10° lower | explained by a 15–25 nm dielectric (oxide) film |
| Diattenuation D | 3 to almost 30× larger; about 0.1 already at 20°, where Al gives < 0.01 | at 20°–30° within the instrument error seen on black glass; at 60°–79° (0.22–0.30 too high) beyond it, but the metals' diattenuation disagrees with Fresnel in general ([QUIRKS.md](QUIRKS.md)): not a property of the sample |
| Depolarization Δ | 0–0.17 on the ridge, above 0.5 in the tails | multiple scattering from roughness; 0 for ideal Fresnel |

What the 32 PCA components could not say, "where does the sample deviate from the ideal", the decomposition separates:
R is the film, D is anisotropy or noise, Δ is roughness.
