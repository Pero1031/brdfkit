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

## Caveats

- MERL's RGB is a camera-calibrated colour (channel scales 1 : 1.15 : 1.66); the hpBRDF RGB here is sRGB with an equal-energy white. Hue components (PC4 and above) therefore
  compare two different colour conventions; the aluminium is nearly grey in the specular region, but the off-peak bins are bluish (median R/B ≈ 0.56 over the φ_d = 90° slice; cause not checked).
- One hpBRDF material only; the other 12 have not been projected. A single point does not say how well the MERL basis generalizes to measured
  hpBRDF metals in general.
- The MERL basis comes from 100 materials of which about a dozen are metals; the metal subset PCA (next item of the roadmap) is the natural follow-up.
- Distances are in feature-norm units of the log (or nielsen) mapping on the common mask; they are comparable between materials of one model, not between models.
