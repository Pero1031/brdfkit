# Editing a measured hpBRDF in the MERL latent space

Goal: materials that were never measured, made by moving a measured hpBRDF through a latent space, while keeping what
only hpBRDF has: 68 spectral bands and full Mueller matrices. Implementation in `brdfkit/analysis/latent_edit.py`;
`examples/latent_demo.py` runs the whole chain from the downloaded data ([WALKTHROUGH.md](WALKTHROUGH.md)).

```bash
python -m brdfkit edit out/factor/Aluminum_s0.npz --model out/merl_analysis/log/merl_pca.npz --out out/latent/al_pc1_-1.npz --pc 1:-1 --keep-albedo 45
python -m brdfkit compress out/latent/al_pc1_-1.npz out/latent/al_pc1_-1_pca.npz --step 4
python -m brdfkit export   out/latent/al_pc1_-1_pca.npz out/latent/al_pc1_-1.hpbpca
```

## Why it works this way

- **hpBRDF alone is too small for a latent space**: 14 materials. MERL has 100, but only RGB and no polarization.
- **The factor model splits the two** ([FACTOR.md](FACTOR.md)): m00 = A(angles) × S(θ_d, λ), m_ij/m00 = N(θ_d, λ).
  A is a grey lobe in luminance units, the same kind of object as a MERL BRDF, so it can live in the MERL latent space;
  S and N stay as measured.
- **Only the change is applied.** The MERL basis does not represent the hpBRDF lobe exactly (projection RMSE 0.41 in log
  units for aluminium, and its flat specular core is a 7.6σ outlier on PC7, [HPBRDF_ON_MERL.md](HPBRDF_ON_MERL.md)).
  Replacing A by its reconstruction would turn it into a MERL-looking lobe. Instead, with z the coefficients of A and dz
  the edit, A' = A · exp(mean over RGB of basis · dz): the decoded difference is applied in the log domain of the 'log'
  mapping, and everything the basis does not capture is kept.

## Steps

1. A is put on the strided MERL grid: φ_d is folded to [0°, 180°) by averaging the four equivalent angles
   (φ, −φ, φ − 180°, 180° − φ), and θ_d and θ_h are taken at the nearest bin (the hpBRDF grid contains the MERL grid).
2. z = masked least-squares coefficients on the first k = 16 components.
3. dz: `--pc P:T` adds T·σ_P along component P, where σ_P is the spread (standard deviation) of the 100 MERL materials
   along it; `--toward NAME:T` moves the fraction T of the way to a MERL material.
4. The log change is interpolated back to the hpBRDF grid (φ_d and θ_d linearly, θ_h linearly in √θ_h, the coordinate in
   which both grids are uniform).
5. `--keep-albedo θ` rescales A' so that its directional albedo at incidence θ equals the original's: only the shape
   changes. (S has luminance 1, so the luminance albedo is the cosine integral of A, estimated by Monte Carlo.)

## Aluminium on the MERL components

z of the aluminium lobe (floor removed, in units of the MERL spread): PC1 −0.5, PC2 −1.1, PC3 −1.5, …, PC7 −7.6.
The components, as read from the MERL analysis ([MERL_PCA.md](MERL_PCA.md)): PC1 specular versus diffuse, PC2 brightness,
PC3 lobe width. Shifts tried, with the albedo at 45° kept (lobe half width: θ_h where A falls to half its peak at
φ_d = 90°, θ_d = 30°; aluminium itself: 1.9°, peak 18):

| Shift | Half width | Peak | Albedo 0° / 75° | Nearest MERL materials after the shift |
|---|---|---|---|---|
| PC1 −1 | 0.28° | 182 | 0.33 / 0.06 | aluminium, nickel, blue-metallic-paint2 |
| PC1 −0.5 | 0.71° | 63 | 0.36 / 0.07 | nickel, aluminium, gold-metallic-paint3 |
| (aluminium) | 1.9° | 18 | 0.37 / 0.09 | |
| PC1 +0.5 | 3.6° | 4.4 | 0.37 / 0.13 | silver-metallic-paint, two-layer-gold, silver-paint |
| PC1 +1 | 6.4° | 0.9 | 0.36 / 0.17 | silver-paint, silver-metallic-paint2, pearl-paint |
| PC3 −1.5 | 4.4° | 3.1 | 0.36 / 0.08 | nickel, two-layer-gold, gold-metallic-paint3 |
| PC3 +1.5 | 0.28° | 84 | 0.35 / 0.15 | violet-acrylic, pink-jasper, two-layer-gold |

Rendered in GoLD_rayt (spectral build, environment light, five spheres per row, the unedited aluminium in the middle):
along **PC1** the sphere goes from a sharp mirror-like metal to a satin, paint-like finish, clearly; along **PC3** the
highlights grow or shrink, more subtly. The degree of linear polarization keeps the aluminium's rim pattern in every
variant, since N is unchanged. These are materials with the colour and polarization of oxidized aluminium and a lobe that
the dataset does not contain.

## Limits

- **Large shifts leave the plausible materials.** ±2σ along PC1 or PC3 without rescaling gives albedos above 1 and peaks
  300 times higher; within ±1–1.5σ and with `--keep-albedo` the results stay below 1. `edit` prints the albedo at 0°, 45°
  and 75° and warns above 1.
- **S and N are fixed.** A rougher surface depolarizes more in reality (multiple scattering); the edit does not model that.
- **A diagonal streak** from the highlight shows on the broadened variants (PC1 +0.5, +1, PC3 +1.5) and faintly on the
  original. Whether it is an amplified φ_d asymmetry of the measurement or a side effect of the interpolation is not known.
- The hpBRDF peak is flattened by the measurement's angular resolution ([FLOOR.md](FLOOR.md), [HPBRDF_ON_MERL.md](HPBRDF_ON_MERL.md));
  sharpening edits sharpen that flattened peak, they do not recover detail that was not measured.
