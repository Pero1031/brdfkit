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
R is smooth along both wavelength and angle and is the most reliable quantity in this analysis.

## Result 2: the diattenuation D is explained by neither Fresnel nor the film

| θ | measured D, 550 nm | bare Al | 20 nm film | measured D, 800 nm | bare Al |
|---|---|---|---|---|---|
| 10° | 0.099 | 0.001 | 0.001 | 0.193 | 0.002 |
| 45° | 0.187 | 0.031 | 0.031 | 0.252 | 0.050 |
| 70° | 0.321 | 0.099 | 0.073 | 0.432 | 0.167 |

The measured D is 3–10 times the theory and stays at 0.1–0.2 even near normal incidence.
On an isotropic surface s and p are indistinguishable at normal incidence, so D must be 0 there; the residual is not material
physics but indicates either

- an anisotropic sample (rolled or brushed aluminium reflects differently along and across the grooves, giving diattenuation even at normal incidence), or
- a systematic error of the polarimeter (averaged out by the dataset's isotropy assumption).

D is also noisy along wavelength (`lu_chipman_spectral.png`, left).
The film barely changes D (the dashed film curve is almost identical to the dotted bare one), so **the thickness information comes from R alone**.

## Result 3: the depolarization Δ is small and grows away from the specular ridge

- In the specular configuration Δ = 0.00–0.17, DI = 0.83–1.00: a reasonable range for a metal of low roughness.
- Δ grows with angle (550 nm: 0.05 at 10°, 0.17 at 70°).
- Off the specular ridge Δ rises to 0.25–0.5 and exceeds 0.7 where the light is near the normal and the viewer near grazing
  (`lu_chipman_map550.png`, right): multiple scattering and the tail component destroy polarization.
- Around 800 nm Δ goes negative (−0.04) and DI exceeds 1. The matrix slightly leaves the physically admissible range
  (sum of squares of |m_ij| > 3 m00²), presumably from measurement noise plus compression error.
  This argues for a physically constrained compression (projection onto passive matrices).

## Note: symmetry of the table

At θ_h = 0, isotropy means the Mueller matrix should not depend on φ_d, yet the φ_d = 0 and φ_d = π slices differ by up to
11% at 20° and 1% at 60°. Off the specular ridge, R is also asymmetric in θ_i and θ_o. Both show that the table is not
strictly isotropic or reciprocal; harmless for rendering, but analyses should average over them.

## Summary

| Quantity | Measured vs ideal Al | Interpretation |
|---|---|---|
| Retardance R | principal angle 8–10° lower | explained by a 15–25 nm dielectric (oxide) film |
| Diattenuation D | 3–10× larger, non-zero at normal incidence | sample anisotropy or a measurement systematic; not material physics |
| Depolarization Δ | 0–0.17 on the ridge, above 0.5 in the tails | multiple scattering from roughness; 0 for ideal Fresnel |

What the 32 PCA components could not say, "where does the sample deviate from the ideal", the decomposition separates:
R is the film, D is anisotropy or noise, Δ is roughness.
