# Editing a measured hpBRDF in the MERL latent space

Goal: materials that were never measured, made by moving a measured hpBRDF through a latent space, while keeping what
only hpBRDF has: 68 spectral bands and full Mueller matrices. Implementation in `brdfkit/analysis/latent_edit.py`;
`examples/latent_demo.py` runs the whole chain from the downloaded data ([WALKTHROUGH.md](WALKTHROUGH.md)).

```bash
python -m brdfkit edit out/factor/Aluminum_s0.npz --model out/merl_analysis/log/merl_pca.npz --out out/latent/al_pc1_-1.npz --pc 1:-1 --keep-albedo 45
python -m brdfkit compress out/latent/al_pc1_-1.npz out/latent/al_pc1_-1_pca.npz --step 4
python -m brdfkit export   out/latent/al_pc1_-1_pca.npz out/latent/al_pc1_-1.hpbpca
```

`edit` first symmetrizes the lobe over φ_d (`--sym median`, the default; see [The streak](#the-streak-a-one-sided-arm-of-the-measured-lobe));
`edit` without `--pc` or `--toward` writes just the symmetrized material, the reference to compare edits with.

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

0. A is made equal at the four φ_d an isotropic, reciprocal BRDF cannot tell apart (φ, −φ, φ + 180°, 180° − φ):
   `--sym median` (default) takes the mean of the middle two of the four values, `--sym mean` their mean, `--sym none` skips it.
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
PC3 lobe width. Shifts tried, lobe symmetrized (median), albedo at 45° kept (half width: hwhm of the φ_d-averaged
microfacet distribution D at θ_d = 30°, measured from its maximum outward as in [GLOSS.md](GLOSS.md); centre: D at θ_h = 0
over its maximum; peak: A at φ_d = 90°, θ_d = 30°):

| Shift | Half width | Centre | Peak | Albedo 0° / 75° | Nearest MERL materials after the shift |
|---|---|---|---|---|---|
| PC1 −1 | 0.23° | 1 | 181 | 0.31 / 0.06 | aluminium, nickel, blue-metallic-paint2 |
| PC1 −0.5 | 0.60° | 1 | 63 | 0.33 / 0.06 | nickel, aluminium, gold-metallic-paint3 |
| (aluminium) | 1.80° | 1 | 18 | 0.35 / 0.08 | two-layer-gold, gold-metallic-paint3, gold-metallic-paint2 |
| PC1 +0.5 | 2.47° | 0.53 | 4.6 | 0.37 / 0.10 | two-layer-gold, silver-metallic-paint, silver-paint |
| PC1 +1 | 2.98° | 0.21 | 0.9 | 0.39 / 0.14 | silver-paint, silver-metallic-paint2, pearl-paint |
| PC3 −1.5 | 2.55° | 0.26 | 3.1 | 0.32 / 0.06 | nickel, two-layer-gold, gold-metallic-paint3 |
| PC3 +1.5 | 0.20° | 1 | 87 | 0.38 / 0.12 | violet-acrylic, pink-jasper, gold-metallic-paint3 |

(An earlier version of this table gave larger half widths for the widened lobes, 3.6° to 6.4°, from a cruder measure.)
The widening shifts dig a hole in the centre of the lobe: its maximum moves to a ring (θ_h 1.1° at PC1 +0.5, 1.6° at
+1), which a point light shows as a ring-shaped highlight. The narrowing shifts go five to six times below the
instrument's resolution (black glass: 1.15°). For the width alone, `brdfkit gloss` is the better knob ([GLOSS.md](GLOSS.md)).

Without the symmetrization (`--sym none`) the half widths and peaks are the same within 5%, and the albedo at 75° is
0.01–0.03 higher: the arms removed by the median carry a little energy, mostly toward grazing.

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
- **Symmetric φ_d structure stays.** The symmetrization removes one-sided structure only (next section).
- **Two brightness scales.** Neither MERL nor the hpBRDF tables as delivered are on an absolute scale, and the two differ
  ([CALIBRATION.md](CALIBRATION.md)). The coefficient on PC2 (brightness), the nearest materials and `--toward` compare
  absolute levels; `--pc` with `--keep-albedo` does not.
- The hpBRDF peak is flattened by the measurement's angular resolution ([FLOOR.md](FLOOR.md), [HPBRDF_ON_MERL.md](HPBRDF_ON_MERL.md));
  sharpening edits sharpen that flattened peak, they do not recover detail that was not measured. Widening edits turn the
  flat peak into a ring (above). The lobes are the instrument's response convolved with the surface ([GLOSS.md](GLOSS.md)).

## The streak: a one-sided arm of the measured lobe

Broadened edits (PC1 +0.5, +1, PC3 +1.5) showed a faint diagonal streak from the highlight, and the unedited aluminium
showed it faintly too. It is in the measured tables, not made by the edit:

- Under one distant light and a distant camera, θ_d is constant over a sphere, so the highlight is a slice A(φ_d, θ_h) at
  that θ_d. At θ_d = 25° the aluminium slice has an arm at φ_d ≈ −22° that reaches θ_h = 20°. It is 3.5 times brighter
  than its mirror image at +22°, where an isotropic BRDF has the same value. A one-light render shows it as a streak off the
  highlight's mirror-symmetry line. Broadening edits lift the shoulders (θ_h 5°–20°) where the arm lies, so it becomes visible.
- All four downloaded hpBRDF materials have such arms, at material-specific angles. In the table, the lobe is the CIE Y of
  m00, averaged over θ_h 6°–15°, with the black-glass floor removed; "arm" is its maximum over φ_d:

| Material | Arm φ_d (θ_d 21° → 45°) | Arm / its mirror, θ_d 25° | Arm / median over φ_d, θ_d 25°: measured → symmetrized (median) |
|---|---|---|---|
| Aluminum | −28° → 0° (turns with θ_d) | 3.5 | 2.6 → 1.1 |
| SUJ2 | −55° → −48° (nearly fixed) | 8.0 | 6.6 → 1.3 |
| Fake gold | −150° and +35° (a φ, φ + 180° pair) | 4.6 | 2.7 → 1.4 |
| Black glass | −120° and +64° (a thin line) | 47 | 46 → 12 |
| MERL aluminium, chrome, black-obsidian | — | 1.0–1.7 | 1.2–1.7 (MERL is stored folded) |

- Black glass is exaggerated by the floor removal: its shoulder is close to the floor, so the median left after the
  subtraction is tiny. Before the removal, its arm is 2.2 times the median and 0.4% of the peak.
- The direction differs between materials while the instrument geometry is the same, so the arms are not one fixed
  instrument artifact (such as a ghost at a fixed place). Per-sample anisotropy (polishing or machining marks) or
  per-sample optics are possible. SUJ2's pattern, one diagonal bright and the perpendicular one dark, is what an
  anisotropic surface gives. The white smooth plastic, a moulded surface, has no arm (1.1 times its mirror image) and the red rough
  plastic a weak one (1.4) ([QUIRKS.md](QUIRKS.md)), so the arms are most likely the samples' own.
- Separately, at θ_d ≥ 65° all four materials are 3–7 times brighter at φ_d = 0° than at 180° (light and camera swapped).
  This reciprocity break is common to all materials, so it more likely comes from the instrument or the processing.

The median of the four equivalent values removes a one-sided arm: the aluminium lobe ends up flatter in φ_d than MERL
aluminium (arm / median 1.1–1.5 against 1.6–2.0 at θ_d 22°–30°). Structure present at both φ and φ + 180° is two of the
four values, so it is halved and copied to −φ and 180° − φ instead of removed: Fake gold and SUJ2 keep a weaker symmetric
pattern, and black glass's line becomes an X. `--sym mean` spreads every arm over the four angles instead (aluminium
1.2–1.6). Rendered in GoLD_rayt with one small light, the streak is gone from the symmetrized aluminium and from its PC1 +1 edit.
A faint, symmetric star of short rays around the PC1 +1 highlight remains (the halved φ + 180° partner of the arm). With the
environment light, the tail under the highlights of PC1 +0.5 and +1 is gone.
