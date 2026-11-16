# PCA compression (`brdfkit compress` / `eval` / `export`)

Compresses a measured BRDF table into a short latent vector per angular bin plus a small basis.
The 13 GB hpBRDF never has to fit in memory: it is streamed one φ_d slice (about 36 MB) at a time.

```bash
python -m brdfkit info     data/kaist/hpbrdf/Aluminum.hpbrdf
python -m brdfkit compress data/kaist/hpbrdf/Aluminum.hpbrdf out/pca/Aluminum_pca.npz --k1 8 --k2 32
python -m brdfkit eval     data/kaist/hpbrdf/Aluminum.hpbrdf out/pca/Aluminum_pca.npz --sweep
python -m brdfkit export   out/pca/Aluminum_pca.npz out/pca/Aluminum.hpbpca
```

| Option | Meaning | Default |
|---|---|---|
| `--k1` | number of principal components of log m00 along the channel axis (68 wavelengths, or RGB 3) | 8 |
| `--k2` | number of principal components of the normalized Mueller elements m_ij/m00 (channels × 15). 0 for unpolarized tables | 32 |
| `--dtype` | storage type of the coefficients (`float16` / `float32`) | float16 |
| `--step` | φ_d stride during the covariance pass (speed-up) | 1 |

## Representation

Each angular bin is represented by two PCAs.

- **Part A (intensity)** `x1 = log(max(m00, eps))`, one dimension per channel → k1 coefficients.
  The logarithm tames the dynamic range of the specular peak, so the PCA does not spend its components on the peak alone.
- **Part B (polarization, polarized tables only)** `x2 = clip(m_ij / m00, −1, 1)` (ij ≠ 00), channels × 15 dimensions → k2 coefficients.
  A passive Mueller matrix satisfies |m_ij| ≤ m00, so clipping to [−1, 1] is physically justified.

Decoding:

```
m00(ch)   = exp( mean_logm00 + basis_logm00 · c1 )
m_ij(ch)  = m00(ch) · ( mean_ratio + basis_ratio · c2 )[ratio_index]
```

Unmeasured bins (NaN) are excluded from the fit and the evaluation; in the coefficient table they are filled with the
coefficients of the nearest measured bin of the same φ_d slice. The `valid` mask stays in the table, so a renderer can
substitute its own fill.

## Output

- `*.npz`: coefficients `coef_logm00` `coef_ratio`, `valid`, means and bases, angular axes, `channels`, `channel_kind`, `polarized`, `meta` (JSON).
  The reference decoder is `brdfkit/decode.py` (`PCADecoder`).
- `*_basis_full.npz`: every eigenvector and the explained-variance ratios (for `eval --sweep`).
- `*.hpbpca` (`export`): a flat binary for renderers (format v2; layout at the top of `cpp/hpbrdf_pca.hpp` and in `brdfkit/export/hpbpca.py`).
  The k1 + k2 coefficients of each bin are stored contiguously; `cpp/hpbrdf_pca.hpp` does the reading and the decoding.

## Results

**hpBRDF Aluminum (k1 = 8, k2 = 32, float16)**

| Item | Value |
|---|---|
| Source file | 13.0 GB |
| Compressed `.npz` / `.hpbpca` | 231 MB |
| Time | covariance 209 s + encoding 105 s |
| Cumulative explained variance | log m00: 99.85% with 8 components; ratios: 99.4% with 32 |
| Relative RMSE of m00 per wavelength | 10.3% (essentially measurement noise: the residual is white along wavelength and matches the noise level estimated from second differences) |
| Relative error of m00 integrated over the 68 bands | median 0.2–0.4%, 99th percentile about 5% |
| RMSE of the normalized Mueller elements | 0.007–0.023 |
| Directional albedo at normal incidence | source table 0.353, decoded 0.351 |

**MERL aluminium (RGB, k1 = 3, float16)**: 9.7 MB, relative RMSE of m00 0.23% (3 channels held by 3 components, so only the float16 rounding remains).

**Known limitation**: right next to the specular peak the decoded value can be 10–18% below the source table
(interpolation in coefficient space = interpolation in log space). If the peak matters, raise the θ_h resolution
or look up the source table around the peak.
