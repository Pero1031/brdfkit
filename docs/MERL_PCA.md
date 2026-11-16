# Inter-material PCA of the 100 MERL materials (`brdfkit merl-pca`)

Reproduces, on the 100 materials of the MERL BRDF database, the analysis that opened the data-driven reflectance model of
Matusik et al. (2003). Each material's whole BRDF table becomes one feature vector, the materials are the samples of a PCA,
and the output is the common basis, the distribution of the materials in it, and the generalization error when a material
is held out of the fit (leave-one-out). Implementation in `brdfkit/analysis/merl_pca.py`.

```bash
python -m brdfkit merl-pca data/merl --out out/merl_analysis/log
python -m brdfkit merl-pca data/merl --out out/merl_analysis/nielsen --mapping nielsen
```

| Option | Meaning | Default |
|---|---|---|
| `--stride` | grid subsampling, one integer or `phi_d,theta_d,theta_h` | 2 (90 × 45 × 45 bins, 1.7 MB per material) |
| `--k` | number of basis vectors stored (coefficients and explained variance are kept for every component) | 16 |
| `--mapping` | feature mapping. `log` = log(max(f, ε)) (Matusik 2003); `nielsen` = log((f cos θi cos θo + ε) / (median over materials + ε)) (Nielsen, Jensen, Ramamoorthi 2015) | log |
| `--eps` | ε | log: 1e-5, nielsen: 1e-3 |
| `--no-plot` | skip the figures | |

## Method

- The bins measured in every material (in MERL the 79% that the geometry allows, common to all materials) × 3 RGB channels
  are concatenated: one material = one vector (D ≈ 440,000 at stride 2).
- No D × D covariance is formed. The Gram matrix (100 × 100) of the centred vectors is eigendecomposed; every eigenvalue
  (= the explained variance of every component), every coefficient of every material, and each material's residual after
  truncation to k components follow exactly from that 100 × 100 matrix.
- Leave-one-out: a PCA on the other 99 materials, then the residual of material i projected onto that basis. This too is
  pure Gram-matrix algebra (100 eigendecompositions of 99 × 99 matrices). The gap to the reconstruction error is how much
  the basis matched a material only because it had seen it.
- Only the first k basis vectors are expanded back to D dimensions and stored. `MaterialPCA.project(values)` maps another
  table on the same grid to coefficients; `reconstruct(coef)` turns coefficients back into linear BRDF values on the strided grid.

## Output (under `--out`)

- `merl_pca.npz`: `names`, `mean [D]`, `basis [D, k]`, `coef [n, n]` (all components), `eigvals`, `evr`, `valid [P', T', H']`, the strided axes,
  `mapping`, `eps`, `ref` (nielsen medians), `cosw` (cos θi cos θo), `recon_rmse [n, n+1]`, `loo_rmse [n, n−1]`, `meta`.
  Load with `MaterialPCA.load(path)`.
- `merl_pca_report.txt`: cumulative explained variance and errors per k; per material its coefficients, errors and nearest neighbour.
- `merl_pca_variance.png` (explained variance and errors vs k), `merl_pca_scatter.png` (materials in PC1–PC2 and PC3–PC4; colour = specular peak / median),
  `merl_pca_basis.png` (mean and PC1–6 as θh × θd images at φd = 90°), `merl_pca_recon.png` (reconstruction of the material with the lowest LOO error and the three with the highest).

## Results (100 materials, stride 2, 2026-11-16)

| k | cum. variance log / nielsen | log-RMSE reconstruction, median | LOO median / p90 / max | nielsen: reconstruction / LOO median |
|---|---|---|---|---|
| 1 | 0.588 / 0.611 | 1.01 | 1.03 / 1.51 / 2.68 | 0.84 / 0.86 |
| 2 | 0.814 / 0.814 | 0.64 | 0.66 / 1.15 / 1.35 | 0.55 / 0.57 |
| 3 | 0.893 / 0.894 | 0.49 | 0.51 / 0.82 / 1.19 | 0.44 / 0.46 |
| 5 | 0.965 / 0.964 | 0.30 | 0.32 / 0.47 / 0.67 | 0.25 / 0.28 |
| 8 | 0.980 / 0.981 | 0.22 | 0.24 / 0.38 / 0.57 | 0.18 / 0.20 |
| 16 | 0.993 / 0.993 | 0.13 | 0.17 / 0.28 / 0.46 | 0.11 / 0.13 |
| 32 | 0.997 / 0.998 | 0.09 | 0.13 / 0.25 / 0.42 | 0.07 / 0.10 |
| 64 | 0.9995 / 0.9996 | 0.04 | 0.10 / 0.21 / 0.39 | 0.03 / 0.08 |

The log-RMSE is in natural-log units: 0.1 means roughly 10% error per bin. Timing: loading 41 s, PCA including LOO 97 s.

**Explained variance versus error.** Five components explain 96% of the variance, yet the per-bin log error is still 30% (median).
The variance is dominated by the big contrasts "metal vs diffuse" (PC1, 59%) and "dark vs bright" (PC2, 22%); the fine shape of
the specular lobe lives in the remaining few percent. Rendering-grade accuracy (about 0.1 in log) needs 16–32 components.

**Generalization.** With the material held out, k = 16 gives a median of 0.17 and a worst case of 0.46 (grease-covered-steel).
The gap to the reconstruction error is 0.03 at the median but 0.1–0.3 for metals (steel 0.22 → 0.40, chrome-steel 0.27 → 0.43).
Only a dozen of the 100 materials are metals, so the basis does not carry enough variety of specular-peak shapes.
Worst LOO: grease-covered-steel, chrome-steel, steel, black-obsidian, green-acrylic, chrome, aluminium, ss440, tungsten-carbide.
Best: pink-fabric (0.07), gold-paint, white-fabric, neoprene-rubber, pink-felt, polyurethane-foam. Fabrics and paints are well represented.

**Meaning of the components** (`merl_pca_basis.png`, `merl_pca_scatter.png`).

- PC1 (59%): specular vs diffuse. Lowers θh < 10° and raises everything else. Metals are negative (chrome −1900, tungsten-carbide, steel, aluminium −1170),
  fabrics positive (pink-fabric, white-fabric, polyethylene +1190).
- PC2 (22%): overall brightness and grazing behaviour. black-obsidian, black-fabric and black-phenolic positive; alumina-oxide, white-acrylic and white-marble negative.
- PC3 (8%): lobe width. Lowers θh 10–30° and raises θh above 50°. specular-*-phenolic and maroon-plastic positive, metallic paints and nickel negative.
- PC4 (5%): hue. blue-acrylic, green-plastic and specular-blue-phenolic positive; red-metallic-paint, specular-red/orange-phenolic negative
  (invisible in the channel-averaged images).
- PC5 onwards: 2% or less; PC6 is 0.6%.

**Nearest neighbours.** In the 16-d coefficient space the nearest neighbours are the name pairs (red-fabric ↔ red-fabric2, silver-paint ↔ pearl-paint,
cherry-235 ↔ ipswich-pine-221, gold-metallic-paint2 ↔ two-layer-gold, chrome ↔ tungsten-carbide). aluminium's neighbour is brass; black-obsidian's is specular-black-phenolic.

**The nielsen mapping.** Multiplying by cos θi cos θo and dividing by the per-bin median over materials lowers the error by 15–25% at every k
(k = 16: reconstruction 0.11, LOO 0.13). The cosine tames the large grazing values, the median removes the shape common to all materials,
and the PCA only has to explain the differences between materials. The price is losing the boundary bins where cos θi cos θo = 0 (common bins 77%).

**Error in linear space.** Reconstructing a fitted material with k components gives a relative RMSE of f cos θi cos θo of 0.17 (median) and 0.32 (p90) at k = 16.
The immediate neighbourhood of the specular peak dominates, so this is a different purpose and a different accuracy from the per-bin PCA of
`docs/PCA.md` (relative RMSE of m00 0.2%). This tool looks at the space of materials; to hand one material to a renderer, use the per-bin PCA.

## Connecting to hpBRDF

`MaterialPCA.project` takes a table on the MERL grid (180 × 90 × 90, RGB). Placing hpBRDF Aluminum (68 wavelengths, Mueller matrices,
91 × 91 × 361 grid) into this space needs a conversion that reduces m00 to RGB and resamples onto the MERL grid. Not implemented yet.
