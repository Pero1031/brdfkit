# The instrument floor of the hpBRDF tables

Far from the specular lobe the hpBRDF tables never drop below about 0.006 sr⁻¹ at 550 nm, even for smooth black glass,
whose true reflectance there is close to zero (MERL black-obsidian: 6e-5). This document shows that this floor is an
additive, unpolarized, per-wavelength offset of the instrument, and describes the tools that estimate and remove it.
Implementation in `brdfkit/analysis/floor.py`.

```bash
python -m brdfkit floor    data/kaist/hpbrdf/Black_glass.hpbrdf --out out/floor/Black_glass_floor.npz
python -m brdfkit compress data/kaist/hpbrdf/Aluminum.hpbrdf out/pca_floor/Aluminum_floor_pca.npz --floor out/floor/Black_glass_floor.npz
python -m brdfkit merl-project data/kaist/hpbrdf/Black_glass.hpbrdf --model out/merl_analysis/log/merl_pca.npz --out out/hpbrdf_on_merl_Black_glass --floor out/floor/Black_glass_floor.npz
```

## Region

"Far from the lobe" is θ_h ≥ 40° and θ_d = 25°–75°. The θ_d limits keep the region inside the measured range: the light
source of the hpBRDF setup is rotated through 40°–160° only, so θ_d below about 20° and above about 80° is filled by Gaussian
inpainting (Moon et al. 2025, Sec. 4; see [LU_CHIPMAN.md](LU_CHIPMAN.md)). Every 12th φ_d slice is read (31 of 361), which
gives 17,421 bins per material and takes about 20 s per material.

## Evidence

**Same spectrum in every material, unrelated to the material's colour.** Median m00 in the region, relative to 550 nm:

| | 414 nm | 446 | 502 | 550 | 598 | 646 | 702 | 798 | 902 | 950 |
|---|---|---|---|---|---|---|---|---|---|---|
| Black glass | 4.97 | 3.00 | 1.33 | 1 | 0.73 | 0.51 | 0.52 | 0.78 | 1.30 | 1.98 |
| Aluminum | 3.15 | 2.09 | 1.23 | 1 | 0.85 | 0.72 | 0.70 | 0.80 | 1.12 | 1.51 |
| Fake gold | 4.38 | 2.66 | 1.28 | 1 | 0.78 | 0.58 | 0.61 | 0.87 | 1.35 | 1.99 |
| SUJ2 | 4.20 | 2.61 | 1.28 | 1 | 0.79 | 0.60 | 0.61 | 0.82 | 1.28 | 1.84 |

Fake gold reflects little blue (its specular peak at 414 nm is 0.27 of the 550 nm value), yet its floor at 414 nm is 4.4
times the 550 nm value. So the floor is not flare or ghosting of the material's own highlight. It is largest where the
lamp × sensor sensitivity is lowest (blue and NIR), which is what a constant added to the raw intensities looks like after
the radiometric calibration divides by that sensitivity.

**No polarization structure.** On black glass every median |m_ij / m00| in the region is below 0.036 (0.021 at 550 nm).
A constant added to every raw intensity, independent of the two QWP angles, is fitted by the Mueller reconstruction
(Eq. 5 of the paper) as m00 alone. Light from the lamp scattered in the room would carry the illumination's polarization
modulation and show up in other elements, so stray lamp light is unlikely; a dark/bias residual or ambient light is not.

**Angular shape.** On black glass (550 nm) the floor does not depend on the cosine of one direction (0.0059–0.0066 over all
bins of cos 0.02–1) and grows about 2× towards grazing in the other (0.0042 at cos 0.9–1, 0.0086 at cos 0.02–0.1).
It is not 1 / cos θ_i, which a constant radiance offset converted to BRDF units would give. A constant per wavelength is
the best simple model (robust spread of log m00 0.44, against 0.9–1.3 for the cosine-weighted variants).

**Metals contain real scattering as well.** The floor of aluminium (0.0117 at 550 nm) has a reflection-like diagonal
(0.18, −0.15, −0.16). Taking the black-glass floor (0.0064) as the instrument offset, aluminium's own far-field scattering is
about 0.0055 and keeps about 35% of the polarization. Fake gold and SUJ2 (0.0075–0.0078, diagonal 0.05–0.08) have little.

## Removal

The black-glass floor is taken as an upper bound of the instrument offset c(λ) (median over the region; `--stat` chooses
another percentile) and removed from m00 only:

m00' = max(m00 − c(λ), keep · m00),   keep = 0.02 by default.

Since the offset is unpolarized, the exact correction leaves the other 15 elements alone, which raises every m_ij / m00 by
m00 / m00'. Where the real signal is much weaker than the offset (black glass far from the lobe) this mostly amplifies noise,
and the compression then spends its polarization components on that noise. So the other elements are scaled by
min(1, max_gain · m00' / m00): no normalized element grows by more than `max_gain` (4 by default, `--floor-max-gain 0` for no
cap). Finally m00' is raised to at least |m_ij| and the norms of the first row and column, so the diattenuation and
polarizance stay ≤ 1. In the lobe (m00 ≫ c) nothing but m00 changes.
`compress --floor` and `eval --floor` apply this per wavelength while streaming the table. `merl-project --floor` applies
the same subtraction to the cached MERL-grid RGB table (the offset is taken to RGB with the same band weights; identical to
subtracting per wavelength except where the clamp acts) and writes its reports to `<out>/<mapping>_floor`.

## Effect on the projection onto the MERL PCA

Offset in RGB: (0.0040, 0.0058, 0.0193). Floor (median far from the lobe, green) before → after: aluminium 0.0138 → 0.0080,
black glass 0.0064 → 0.0005, Fake gold 0.0078 → 0.0020, SUJ2 0.0077 → 0.0019.

| | log: nearest (before → after) | log residual k=16 | nielsen residual k=16 | z of PC7 (log) |
|---|---|---|---|---|
| Black glass | blue-acrylic, green-acrylic, … (black-obsidian 100th) → **black-obsidian, specular-black-phenolic, black-phenolic** | 0.38 → 0.84 | 0.29 → 0.41 | −6.1 → −10.7 |
| Aluminum | two-layer-gold, two-layer-silver, … → two-layer-gold, gold-metallic-paint3, gold-metallic-paint2, nickel | 0.32 → 0.50 | 0.28 → 0.34 | −7.0 → −7.6 |
| Fake gold | green-acrylic, violet-acrylic, … → silicon-nitrade, red-metallic-paint, ss440 | 0.51 → 0.90 | 0.37 → 0.41 | −8.0 → −11.2 |
| SUJ2 | two-layer-gold, gold-metallic-paint2, … → nickel, color-changing-paint2, grease-covered-steel | 0.37 → 0.81 | 0.30 → 0.41 | −6.5 → −7.1 |

- The neighbours become plausible: black glass now sits next to black-obsidian and the black phenolics; steel (SUJ2) next to
  nickel and grease-covered steel.
- The residuals grow, most in the log mapping. What is left after the subtraction in the far field is small and noisy
  (and clamped at 2% of the measurement), and the log mapping weighs those small values heavily. The nielsen mapping, which
  compresses the low end, changes less.
- The PC7 outlier (the flat specular core, θ_h < 3°) does not go away and gets stronger in log: it comes from the angular
  resolution of the measurement, not from the floor.

## Effect on the compression

`compress --k1 8 --k2 32`, median RMSE of the normalized Mueller elements m_ij / m00 after decoding, by distance from the
specular ridge (12 φ_d slices; reference = the table with the same floor treatment):

| | θ_h < 5° | 5°–20° | ≥ 40° | ratio variance in 8 components |
|---|---|---|---|---|
| Black glass, no floor removal | 0.010 | 0.014 | 0.014 | 0.98 |
| Black glass, floor removed, no cap | 0.021 | 0.058 | 0.100 | 0.71 |
| Black glass, floor removed, cap 4 | 0.012 | 0.034 | 0.043 | 0.96 |
| Aluminum, no floor removal | 0.008 | 0.010 | 0.013 | 0.98 |
| Aluminum, floor removed, no cap | 0.014 | 0.012 | 0.030 | 0.93 |
| Aluminum, floor removed, cap 4 | 0.011 | 0.011 | 0.026 | 0.96 |

Without the cap the far-field noise takes over the polarization components and the error near the highlight doubles; with
the cap it is close to the error without floor removal. The m00 error is about the same in all cases (relative RMSE 0.10–0.11).
