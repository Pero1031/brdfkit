# Port bundle for rendering hpBRDF (and the other BRDF tables brdfkit exports) in GoLD_rayt

The KAIST/POSTECH hpBRDF (measured spectral, polarimetric BRDF; 13 GB per material) is PCA-compressed to 231 MB and wrapped
as a GoLD_rayt `Material` with evaluation, sampling and polarization (`scatterMueller`).
It was written against GoLD_rayt's own code, following its API and conventions, and was built, linked and run against
GoLD_rayt's headers and `GGX.cpp` (`test/build_check.bat`, both the RGB and the spectral build). Rendering inside GoLD_rayt
itself has not been done here yet. The steps are in [PATCH_NOTES.md](PATCH_NOTES.md).

## Contents

| File | Role |
|---|---|
| `include/Materials/MeasuredHpBRDF.hpp` | GoLD_rayt `Material` implementation (header only): eval / sample / pdf / scatterMueller / aovAlbedo |
| `include/Materials/hpbrdf_pca.hpp` | dependency-free decoder (`namespace hpbrdf`): file reading, Rusinkiewicz coordinates, latent-vector interpolation, decoding |
| `out/pca/Aluminum.hpbpca` (repository-level out/, not tracked) | compressed aluminium (231 MB, latent 8+32 dimensions, float16). Made with `python -m brdfkit export` |
| `scenes/hpbrdf_aluminum.json` | example scene: the measured aluminium sphere next to a Johnson & Christy gold sphere |
| `test/` | a translation unit and batch file that only check that it compiles against GoLD_rayt's headers |

## Mapping GoLD_rayt's conventions

| GoLD_rayt | Decoder (Mitsuba 3 `measured_polarized` conventions) |
|---|---|
| `wo` (toward the viewer), `wi` (toward the light), world space | `to_view`, `to_light`, converted into the local frame `frame::Frame(rec.n)` |
| `eval` returns f_r without the cosine term | the table's m00 also has no cosine term; returned as is |
| Spectrum = 4 wavelengths (spectral build) | `m00(c, λ)` decoded per wavelength |
| Spectrum = RGB (RGB build) | at load time the 68-band m00 is baked into an RGB table through the colour-matching functions (about 3 M bins × 3), normalized to equal-energy white |
| `scatterMueller`: M/M00 in the s basis (s ⟂ plane of incidence, shared by both beams) | rotated from the table's frame (reference vectors in the plane of reflection) with `pol::basisChange`: `basisChange(x_out, s, wo) · M/m00 · basisChange(s, x_in, −wi)` |
| Stokes: (x, y = k × x, k), exp(−iωt) | Mitsuba 3 is also right-handed (x, y, k). Signs verified below |
| Importance sampling | mixture of GGX VNDF (`alphaSample`) 0.9 + cosine 0.1. Unbiased for any α |

## Checking the polarization sign conventions

The table's Mueller matrices were rotated into GoLD's s basis with the formula above and, in the specular configuration (θ_h = 0),
compared with the Mueller matrix that GoLD's own `pol::fresnelAmplitudes` + `setInterface` produce for Al (Rakić n, k).

| Element | Table (450–800 nm, θ = 20–75°) | GoLD's Fresnel | Sign |
|---|---|---|---|
| m01 | +0.06 to +0.48 | +0.005 to +0.22 | agree (s reflects more strongly) |
| m22, m33 | negative (crossing 0 at grazing) | negative | agree |
| m23 | +0.06 to +0.91 | +0.03 to +0.91 | agree |
| m32 | −0.02 to −0.79 | −0.03 to −0.91 | agree |

All signs agree, so no sign flip of S3 (circular polarization) is applied.
The differences in magnitude come from the measured aluminium not being an ideal Fresnel surface (oxide film and roughness make
m01 larger and |m22| smaller), which is exactly the point of using measured BRDFs.

## Checks (`test/check_rgb.exe`, `test/check_spectral.exe`, run from the bundle root)

| Item | RGB build | Spectral build (450/550/650/800 nm) |
|---|---|---|
| Loading (including the RGB bake) | 0.85 s | 0.35 s |
| f_r in the specular direction (θ = 45°) | 11.8 / 13.7 / 12.7 | 12.2 / 13.7 / 11.6 / 11.3 |
| f_r off-specular | 0.018 / 0.021 / 0.034 | 0.033 / 0.021 / 0.019 / 0.016 |
| `scatterMueller` (specular, 450 nm) m01 / m22 / m23 / m32 / m33 | — | +0.127 / −0.792 / +0.428 / −0.396 / −0.874 (matches the Python reference) |
| Sampling E[f cos / pdf] (θ_o = 45°, 20k samples) | 0.231 | 0.267 (450 nm) |
| One eval | 441 ns | 776 ns |

**On the brightness of the aluminium.** The directional albedo of this measured aluminium is 0.35 at normal incidence,
0.26 at 45° and 0.08 at 75°, far below the Fresnel reflectance of ideal Al (about 0.9), and darker at oblique angles.
Nothing was lost in compression (integrating the source table directly gives 0.353, the decoded table 0.351); the dataset's
sample simply reflects that little. It looks dark next to the `conductor` Au, but that is the measurement.

**Reciprocity.** The table is not strictly reciprocal (swapping wo and wi changes off-specular values by 10–20% in places).
Like the Mitsuba plugin, it is used as is. If needed, φ_d and φ_d + π can be averaged at compression time to symmetrize it.

## Performance and memory

| Item | Value |
|---|---|
| Loading | 0.33 s (float16 → float32 expansion), plus a few seconds for the RGB bake in the RGB build |
| RAM | coefficients 456 MB (float32) + RGB table 36 MB |
| One evaluation | about 0.65 µs unpolarized, about 1.0 µs polarized (standalone benchmark; not measured inside GoLD_rayt) |

If memory matters, `hpbrdf_pca.hpp` can be changed to keep the coefficients in float16 and convert during interpolation (halving to about 230 MB).

## Known limitations

- About 25% of the bins, at grazing angles, are unmeasured and hold the value of the nearest measured bin in the same φ_d slice. The `valid` mask stays in the table.
- Right next to the specular peak the decoded value can be 10–18% below the source table (interpolation in coefficient space = interpolation in log space).
- Outside 414–950 nm the edge bands are extended. GoLD_rayt's 360–414 nm falls in this range.
- `scatterMueller` in the RGB build represents R/G/B by 610/550/465 nm.
- The other 13 materials (black glass, rough plastics in several colours, silicone, SUJ2, fake gold, billiard ball) can be fetched
  from Hugging Face and converted to the same format with `python -m brdfkit compress` → `export` (about 5 minutes per material).
