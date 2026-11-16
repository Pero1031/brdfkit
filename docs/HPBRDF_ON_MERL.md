# hpBRDF aluminium in the MERL PCA space (`brdfkit merl-project`)

Projects a measured table onto the inter-material PCA basis of the 100 MERL materials ([MERL_PCA.md](MERL_PCA.md)). It answers
two questions: where does the material sit among the MERL materials, and how well does the MERL basis represent it
(the table is a held-out sample, so the residual is a generalization error). Implementation in
`brdfkit/analysis/merl_project.py`; the colour conversion is in `brdfkit/color.py`.

```bash
python -m brdfkit merl-pca data/merl --out out/merl_analysis/log                       # once: the models (see MERL_PCA.md)
python -m brdfkit merl-pca data/merl --out out/merl_analysis/nielsen --mapping nielsen
python -m brdfkit merl-project data/kaist/hpbrdf/Aluminum.hpbrdf \
    --model out/merl_analysis/log/merl_pca.npz out/merl_analysis/nielsen/merl_pca.npz \
    --out out/hpbrdf_on_merl --reference data/merl/aluminium.binary
```

| Option | Meaning | Default |
|---|---|---|
| `--model` | one or more `merl_pca.npz`; each gets its own subfolder named after its mapping | required |
| `--k` | number of components of the projection (at most the number stored in the model) | 16 |
| `--sym` | `full`: average the φ_d-equivalent slices into each MERL bin; `none`: take the slice at the same φ_d only | full |
| `--reference` | MERL `.binary` of a comparison material (compared with the table, and drawn in the figures) | none |
| `--recompute` | ignore the cached MERL-grid table `<out>/<name>_merl_grid_<sym>.npz` | |
| `--no-plot` | skip the figures | |

## Method

1. **Grid.** The hpBRDF axes contain every MERL grid point to float precision (φ_d every 1° over [−180°, 180°], θ_d every 1°,
   θ_h with the same squared spacing as MERL), so nothing is interpolated: the 180 × 90 × 90 MERL bins are picked by index.
   A table whose axes do not contain the MERL grid is rejected.
2. **Colour.** m00 of the 68 bands (414–950 nm) becomes linear sRGB with the band weights of the renderer
   (`MeasuredHpBRDF.hpp`: Wyman–Sloan–Shirley CIE fit, XYZ → sRGB D65, each channel's weights summing to 1 so that a flat spectrum is
   equal-energy white). A test compares the constants with the header. Negative results are clamped at 0.
3. **φ_d.** For an isotropic, reciprocal, mirror-symmetric surface f(φ) = f(−φ) = f(φ + 180°) = f(180° − φ). With `--sym full`
   the slices of these four angles are averaged into the MERL bin (360 slices feed 180 bins; −180° and +180° count once), and their
   relative spread is recorded. MERL itself stores only φ_d ∈ [0°, 180°).
4. **Projection.** Least squares on the first k basis vectors in feature space (`MaterialPCA.project_masked`). The hpBRDF is unmeasured
   in a thin band of bins at the horizon that MERL measured (3.3% of the log model's features), so the plain full-mask `project` cannot
   be used. The restricted basis stays well conditioned (condition number 1.11).
5. **Evaluation.** Feature residual for k = 1..16; the same residual for the MERL materials held out of the fit (leave-one-out, from the
   model); z-scores of the coefficients (coefficient / its standard deviation over the 100 materials); distances to every MERL material in
   the k-d coefficient space; comparison with the MERL aluminium.

## Output (under `--out`)

- `Aluminum_merl_grid_full.npz`: the table on the MERL grid (`values` [180, 90, 90, 3] linear RGB, `spread`, `meta`). Reused on later runs.
- `<mapping>/report.txt`, `projection.json` (coefficients, residuals, nearest materials).
- `<mapping>/hpbrdf_on_merl_scatter.png` (PC1–PC2 and PC3–PC4 with all MERL materials, the hpBRDF aluminium and the MERL aluminium),
  `..._slices.png` (φ_d = 90° slices of the table, the projection and the MERL aluminium, differences, a θ_h profile),
  `..._error.png` (residual versus k against the MERL leave-one-out band, coefficient z-scores).

## Results (hpBRDF Aluminum, k = 16, stride 2, 2026-11-16)

| | log mapping | nielsen mapping |
|---|---|---|
| feature RMSE (ln units) at k = 1 / 3 / 5 / 8 / 16 | 1.76 / 1.08 / 0.96 / 0.42 / 0.32 | 1.51 / 0.93 / 0.77 / 0.48 / 0.28 |
| MERL materials held out of the fit, k = 16: median / p90 / max | 0.166 / 0.282 / 0.461 | 0.131 / 0.204 / 0.385 |
| MERL aluminium, k = 16: held out / in the fit | 0.313 / 0.244 | 0.249 / 0.182 |
| relative RMSE of f cos θi cos θo at k = 16 | 0.25 | 0.32 |
| nearest MERL materials (distance, RMS per feature) | two-layer-gold (1.24), two-layer-silver (1.30), gold-metallic-paint2 (1.33), gold-metallic-paint3, alum-bronze | two-layer-gold (1.09), two-layer-silver (1.20), gold-metallic-paint3, gold-metallic-paint2, nickel |
| MERL aluminium: distance (RMS) / rank among all materials | 1182 (1.79) / 18th | 990 (1.52) / 15th |
| raw feature RMSE to the MERL aluminium, no PCA | 1.84 | 1.56 |

**The hpBRDF aluminium is not near the MERL aluminium.** At k = 16 the MERL aluminium is the 15th–18th nearest of the 100 MERL materials.
The nearest ones are the metallic paints, two-layer gold and silver, and alum-bronze: rougher and lower-peaked metals.
The two samples differ a lot in the raw tables, before any PCA (feature RMSE 1.6–1.8). At φ_d = 90°, θ_d = 30° (green channel):

| θ_h | 0° | 1.1° | 4.4° | 10° | 28° | 54° |
|---|---|---|---|---|---|---|
| hpBRDF Aluminum | 20.6 | 16.2 | 0.89 | 0.104 | 0.017 | 0.015 |
| MERL aluminium | 438 | 2.8 | 0.108 | 0.016 | 0.0036 | 0.0012 |

The hpBRDF peak is about 20 times lower and flat within the first 0.3°, then the lobe is wider and the floor 5–12 times higher; the median
ratio over all bins is ×3.4 / ×3.6 / ×4.0 in R / G / B. A flat top of that kind looks like the angular resolution of the measurement, but this was not checked against
the measurement description.

**The basis represents it about as well as it represents the MERL aluminium, and worse than the typical MERL material.** The residual at k = 16
(0.32 log, 0.28 nielsen) is the same as the MERL aluminium held out of the fit (0.31, 0.25) and about twice the median over the MERL
materials (0.17, 0.13); it is inside the leave-one-out range of the metals (max 0.46, 0.39). The residual map shows ±0.2 (log10) over most of the hemisphere and
±0.8 at grazing θ_h > 70° and θ_d > 70°. In the log model the residual drops from 0.94 to 0.42 between k = 6 and k = 8:
**component 7 is what this material needs.** Its coefficient is −560, z = −7.1 (the MERL range is −249 to +201), and about 80% of its contribution
comes from θ_h < 3°. That component carries the shape of the specular core, which in no MERL material is as low and flat as in this table. In the nielsen model the same material is
extreme in other components (PC6 z = +4.2, PC10 z = +3.6, PC14 z = −2.9), not in one.

**Symmetry of the measurement.** The relative spread of the four equivalent slices is 0.16 (median) and 0.41 (p90) per bin. Measured on single slices,
the mirror pair f(φ) / f(−φ) differs by 13–19% (median |ln ratio|, no clear sign), the reciprocity pair f(φ) / f(φ − 180°) by 13–14% at φ_d = 60–75° but by 34–49% at φ_d = 20–30°,
where f(φ) is the larger (median +12 to +21%). The MERL aluminium differs by 17–20% between f(φ) and f(180° − φ) with the same measure, so the size is comparable to the
noise of MERL itself; the φ_d-small bias is a property of this table (cause not investigated).

**Sensitivity to `--sym`.** With `--sym none` (no averaging) the coefficients move by 66 of 1145 (log) and 47 of 1084 (nielsen) in norm (about 6% and 4%), the five
nearest materials of the log model are unchanged (the MERL aluminium moves from 18th to 16th), and the residual at k = 16 rises from 0.32 to 0.42 (log) and from 0.28 to 0.36 (nielsen):
the averaging removes noise that the basis cannot represent anyway, and does not move the material.

## Is a metal-only basis better? (`brdfkit merl-subsets`)

The projection above leaves two explanations: the 100-material basis lacks the variety of metal specular shapes, or the table is unlike any MERL
metal. `brdfkit/analysis/merl_subsets.py` refits the PCA on subsets of the MERL materials (log mapping, stride 2) and projects the table on each. The subsets are
11 bare metals and steels (aluminium, alum-bronze, brass, chrome, chrome-steel, nickel, ss440, steel, tungsten-carbide, black-oxidized-steel,
grease-covered-steel), those plus 12 metallic paints and two-layer coatings (23), each with and without the MERL aluminium, and, as the control,
30 random subsets of each size drawn from the other materials. Every residual is measured on the same bins (measured in all 100 materials and in the table).

```bash
python -m brdfkit merl-subsets out/hpbrdf_on_merl/Aluminum_merl_grid_full.npz data/merl --out out/subsets     --all-model out/merl_analysis/log/merl_pca.npz --random 30 [--ablate]
```

ln RMSE of the log BRDF (relative RMSE of f cos θi cos θo in the last column):

| basis | n | k = 5 | k = 8 | k = 9 | largest k (value / linear) |
|---|---|---|---|---|---|
| all materials | 100 | 0.96 | 0.42 | 0.41 | 16: 0.32 / 0.25 |
| metals | 11 | 1.03 | 0.94 | 0.91 | 9: 0.91 / 1.01 |
| metals, no aluminium | 10 | 1.25 | 1.15 | | 8: 1.15 / 0.78 |
| random 11, mean ± std (best of 30) | 11 | 1.01 ± 0.15 | 0.78 ± 0.15 (0.53) | 0.75 ± 0.15 (0.52) | 9: 0.75 / 0.83 |
| metals + finishes | 23 | 1.01 | 0.67 | 0.64 | 16: 0.55 / 0.33 |
| metals + finishes, no aluminium | 22 | 1.04 | 0.73 | 0.70 | 16: 0.60 / 0.37 |
| random 23, mean ± std (best of 30) | 23 | 0.99 ± 0.10 | 0.69 ± 0.12 (0.45) | 0.64 ± 0.11 (0.43) | 16: 0.48 ± 0.06 (0.35) / 0.55 |

- **A metal-only basis does not help.** At k = 9 the 11 metals leave 0.91, against 0.75 ± 0.15 for random 11-material subsets and 0.41 for all 100 materials.
  The 23-material set leaves 0.55 at k = 16, again above its random control (0.48 ± 0.06, 1.3 standard deviations). The explanation "the basis lacks metals" is not supported as stated.
- **Variety, not metalness, lowers the residual.** The residual falls with the number of materials, and the metals are held out better by the full basis too:
  their median leave-one-out error at k = 9 is 0.38 under the 100-material basis and 0.50 under the metal-only basis.
- **The MERL aluminium helps, but only a little.** Removing it from the metal basis raises the residual (0.94 → 1.15 at k = 8), yet the table stays far from the best random subsets (0.53).
- **The drop at k = 7 exists only in the full basis.** The residual falls from 0.94 to 0.42 there; none of the metal sets, and neither mean curve of the random subsets (0.93 → 0.83 and 0.92 → 0.81 from k = 6 to 7), shows a drop of that size.
  The component that carries the specular core is therefore not a property of the metals as a class.
- **No single material carries it.** Refitting the 100-material basis 99 times, each time without one material (`--ablate`), moves the k = 16 residual by at most +0.024
  (median change 0.000, std 0.005; the full basis gives 0.324). The largest increases are nickel (+0.024), tungsten-carbide (+0.013), aluminium (+0.012), brass (+0.010)
  and grease-covered-steel (+0.008), so metals are 5 of the top 6, but each effect is within a few percent of the residual. The core component comes from the
  combined variety of the set, which is why neither 11 metals nor 23 materials reproduce it, and why dropping one material does not remove it.

## Three more hpBRDF materials (`brdfkit merl-compare`)

To see whether the aluminium result is specific to that sample, Fake gold, SUJ2 (bearing steel) and Black glass (13 GB each, from the dataset's Hugging Face
page) went through the same pipeline, each into its own folder, and `brdfkit merl-compare` summarizes the four runs:

```bash
python -m brdfkit merl-project data/kaist/hpbrdf/Black_glass.hpbrdf --model out/merl_analysis/log/merl_pca.npz out/merl_analysis/nielsen/merl_pca.npz     --out out/hpbrdf_on_merl_Black_glass --reference data/merl/black-obsidian.binary
python -m brdfkit merl-compare out/hpbrdf_on_merl out/hpbrdf_on_merl_Black_glass out/hpbrdf_on_merl_Fake_gold out/hpbrdf_on_merl_SUJ2
```

| | Aluminum | Black glass | Fake gold | SUJ2 |
|---|---|---|---|---|
| f at θ_h = 0 (φ_d 90°, θ_d 30°, green) | 20.6 | 2.37 | 28.1 | 10.0 |
| f(θ_h = 0.28°) / f(0) (1 = flat top) | 0.975 | 0.961 | 0.958 | 0.971 |
| floor (median f at θ_h ≥ 40°) | 0.0138 | 0.0064 | 0.0078 | 0.0077 |
| φ_d spread of the equivalent slices (median) | 0.163 | 0.066 | 0.126 | 0.139 |
| residual at k = 16, log / nielsen | 0.32 / 0.28 | 0.38 / 0.29 | 0.51 / 0.37 | 0.37 / 0.30 |
| most extreme component, log (z) | PC7 (−7.0) | PC7 (−6.1) | PC7 (−8.0) | PC7 (−6.5) |
| MERL comparison material, rank by distance (log / nielsen) | aluminium 18 / 15 | black-obsidian 100 / 21 | two-layer-gold 11 / 17 | steel 69 / 45 |

For comparison the MERL materials at the same place: black-obsidian f(0) = 43.7 (falling to 3.0 at 0.28°) and a floor of 1e-4; aluminium f(0) = 438, floor 0.002.

- **The core component is a property of the hpBRDF tables, not of the aluminium.** All four materials lean hardest on PC7 in the log model (z from −6.1 to −8.0), and all four have
  the same flat specular top: f at 0.28° is 96–98% of f at 0°, also for the smooth black glass, whose true peak is far sharper (MERL's obsidian loses 93% over the same range).
  A smooth dielectric cannot be flat-topped by itself, so the flat top is consistent with the angular resolution of the measurement. This is an inference from the four tables; it
  was not checked against the measurement description.
- **A common floor.** Far from the lobe the tables sit at 0.0064–0.0078 for the glass, the gold and the steel, three very different surfaces, while the MERL obsidian is at 1e-4. A glass has no
  diffuse scattering at that level, so a floor of about 0.007 is probably an instrument floor (stray light or noise); only the aluminium is higher (0.014), possibly real scattering.
  This floor is also why hpBRDF tables look about 3 times brighter than MERL in the median ln ratio (×3.1 for the glass, ×3.3 for the steel against steel, ×3.5 for the aluminium):
  most bins are far from the lobe. A constant 0.007 adds about 2% to the albedo (π × 0.007).
- **Dark materials are dominated by the floor in the log mapping.** The black glass is the 100th of 100 by distance to the MERL obsidian in the log model (21st in nielsen), because ln(0.006 / 1e-4) is a
  large constant over most of the hemisphere. Distances of dark materials in the log model say more about the floor than about the lobe shapes.
- **Fake gold is the worst represented** (0.51 log, 0.37 nielsen) and lands next to green and violet acrylic rather than gold paints. Its colour (median ratio to MERL two-layer-gold
  0.39 / 0.49 / 0.73 in R / G / B) shows the colour conventions of the two sets differ, as noted below; hue components of the comparison should not be over-read.
- **Symmetry.** The φ_d spread is lowest for the glass (0.07, the 4-slice average has the least to remove) and 0.13–0.16 for the three metals.

## Caveats

- MERL's RGB is a camera-calibrated colour (channel scales 1 : 1.15 : 1.66); the hpBRDF RGB here is sRGB with an equal-energy white. Hue components (PC4 and above) therefore
  compare two different colour conventions; the aluminium is nearly grey in the specular region, but the off-peak bins are bluish (median R/B ≈ 0.56 over the φ_d = 90° slice; cause not checked).
- Four hpBRDF materials of 14 have been projected, all from one instrument. They agree with each other, which does not separate the instrument from the dataset's processing.
- The subset study uses 11 and 23 materials, so at most 9 and 21 components; the comparison with random subsets of the same size is what makes it fair, but a basis
  of a few dozen materials is small whatever it is made of. Only the log mapping was studied.
- Distances are in feature-norm units of the log (or nielsen) mapping on the common mask; they are comparable between materials of one model, not between models.
