# Walkthrough: from measured hpBRDF tables to edited, renderable materials

This page ties the pieces of brdfkit together: what was done to the hpBRDF data, in which order, why, and how to repeat it
after downloading the datasets yourself. The datasets and everything derived from them stay out of the repository
([LICENSES.md](../LICENSES.md)); `examples/latent_demo.py` rebuilds the demo from the downloads.

## The question

The hpBRDF dataset (Moon et al., SIGGRAPH Asia 2025) gives isotropic BRDFs as 4 × 4 Mueller matrices at 68 wavelengths
(414–950 nm) for 14 materials, 13 GB each. A spectral, polarized renderer can use them directly, but material *editing*
(making a material that was never measured by moving through a latent space) needs many materials, and 14 is too few.
MERL (Matusik et al. 2003) has 100 materials but only RGB and no polarization. The approach here is to split each hpBRDF
into a grey angular lobe, which can be edited in the MERL latent space, and a spectrum and polarization, which are kept
as measured.

## Pipeline

| Step | Command | Document | Result (aluminium unless noted) |
|---|---|---|---|
| 1. Load | `brdfkit info` | [README](../README.md) | every dataset as one `BRDFTable` (Mitsuba `measured_polarized` conventions) |
| 2. Know the measurement | — | [LU_CHIPMAN.md](LU_CHIPMAN.md), [HPBRDF_ON_MERL.md](HPBRDF_ON_MERL.md), [LATENT_EDIT.md](LATENT_EDIT.md#the-streak-a-one-sided-arm-of-the-measured-lobe) | only θ_d ≈ 20°–80° is measured (the rest is inpainted); the specular peak is flat within 0.28° (angular resolution); the lobes have one-sided arms in φ_d that an isotropic material cannot have; every lobe is the instrument's response (black glass: half width 1.1° at θ_d 20°–30°, wider at larger θ_d) convolved with the surface |
| 3. Remove the instrument floor | `brdfkit floor`, `--floor` | [FLOOR.md](FLOOR.md) | an unpolarized offset of about 0.006 sr⁻¹ (550 nm), larger in the blue; without removal the materials render bluish and polarization is diluted |
| 4. Factor the table | `brdfkit factor` | [FACTOR.md](FACTOR.md) | m00 = A(angles) × S(θ_d, λ), m_ij/m00 = N(θ_d, λ): 6.8 MB, 6% m00 and 0.01 polarization error near the specular ridge |
| 4b. Calibrate on black glass | `brdfkit calibrate`, `edit --calib` | [CALIBRATION.md](CALIBRATION.md) | the tables carry 0.47 (θ_d 20°) to 0.08 (75°) of the Fresnel energy, for glass and aluminium alike; dividing by this K brings aluminium within 0–15% of its Fresnel reflectance |
| 5. MERL latent space | `brdfkit merl-pca` | [MERL_PCA.md](MERL_PCA.md) | PCA of the 100 MERL materials (log mapping); PC1 specular/diffuse, PC2 brightness, PC3 lobe width |
| 6. Edit | `brdfkit edit` | [LATENT_EDIT.md](LATENT_EDIT.md) | A symmetrized over φ_d and moved along MERL components, S and N kept, albedo kept: mirror-like to satin aluminium |
| 6b. Swap the optics | `brdfkit swap` | [SWAP.md](SWAP.md) | the lobe shape of one material with the colour, polarization and reflectance of another (aluminium → gold, steel, black glass) |
| 6c. Set the gloss | `brdfkit gloss` | [GLOSS.md](GLOSS.md) | the lobe widened to a target half width or GGX-equivalent roughness by convolving its microfacet slopes (aluminium from hwhm 1.8° to 25°), colour and polarization kept |
| 7. Render | `brdfkit compress`, `export` | [PCA.md](PCA.md), [integrations/gold_rayt](../integrations/gold_rayt/) | `.hpbpca` files for the spectral, polarized renderer GoLD_rayt |

Side studies that shaped these choices: the MERL basis does not lack metals (a metal-only basis is worse, `merl-subsets`);
the four hpBRDF materials checked (aluminium, black glass, Fake gold, SUJ2 steel) share the flat peak and the floor, so
both come from the measurement, not the materials ([HPBRDF_ON_MERL.md](HPBRDF_ON_MERL.md)).

## Repeating it

### Data

| What | Where | Size | Put it in |
|---|---|---|---|
| hpBRDF `Aluminum` (the material to edit) | [Hugging Face](https://huggingface.co/datasets/yunseongmoon/Hyperspectral-Polarimetric-BRDF) (the authors' `download_hpbrdf.py` in [their repository](https://github.com/yunseong0518/hpBRDF) fetches it) | 13.0 GB | `data/kaist/hpbrdf/Aluminum.hpbrdf` |
| hpBRDF `Black glass` (reference for the floor and the calibration; optional but recommended) | same | 13.0 GB | `data/kaist/hpbrdf/Black_glass.hpbrdf` |
| MERL BRDF database, 100 `.binary` files | [Zenodo, doi:10.5281/zenodo.8101680](https://doi.org/10.5281/zenodo.8101680) (CC BY-SA 4.0) | 3.3 GB | `data/merl/` |

File names with spaces are expected with underscores. A symlink or junction into wherever the data lives works.
Each hpBRDF file is 13,010,049,768 bytes; check the size after downloading (an interrupted download is easy to miss).
The hpBRDF release states no license: use it for research and ask the authors before redistributing anything derived from it.

### Run

```bash
pip install -e .[plot]
python examples/latent_demo.py --hpbrdf-dir data/kaist/hpbrdf --merl-dir data/merl --out out/demo
```

It runs steps 3–7 above, skipping any step whose output exists, and writes the GoLD_rayt scenes
`out/demo/scenes/row_pc1_sym-median.json` and `row_pc3_sym-median.json` (five spheres each: the edited aluminium, the
symmetrized original in the middle). `--sym none` skips the φ_d symmetrization (output names then have no `_sym-*` suffix).
`--calib CAL.npz` also applies the black-glass calibration (output names get `_cal`); make the calibration first with
`brdfkit factor` on the black glass and `brdfkit calibrate` ([CALIBRATION.md](CALIBRATION.md)).
With the renderer it also renders them and writes a contact sheet (give `--env`: the scenes then need no file relative to
the renderer's folder):

```bash
python examples/latent_demo.py --renderer <exe of the spectral build> --renderer-cwd <GoLD_rayt>/GoLD_rayt --env <an .hdr environment map>
```

GoLD_rayt is at [github.com/Pero1031/GoLD_ray](https://github.com/Pero1031/GoLD_ray); use its spectral build (the RGB
build bakes the spectrum to RGB at load). `--edits "1:-1,-0.5,0.5,1;3:-1.5,1.5"` chooses the components and shifts
(in units of the spread of the MERL materials), `--material` another hpBRDF material. `--renderer-cwd` is the folder GoLD_rayt
reads its own assets from (it prints a harmless warning about a missing gold table otherwise).

| Step | Time on an 8-thread desktop | Disk |
|---|---|---|
| MERL PCA (reads 100 files) | minutes | 30 MB |
| Floor (black glass, 31 of 361 φ_d slices) | 20 s | — |
| Factor model (two passes over 13 GB) | 7–10 min | 7 MB |
| Each edit | seconds | 7 MB |
| Each bake to `.hpbpca` | 2–3 min | 231 MB |
| Each row of five spheres, 1280 px, 64 spp | 3–4 min | — |

## What is not settled

- The cause of the radiometric error K(θ_d) (an illumination-angle calibration or cosine factors in the tables) and its
  value beyond θ_d 65° ([CALIBRATION.md](CALIBRATION.md)); the spectral errors of the glass do not repeat on aluminium, so S is not corrected.
- The floor's source (dark/bias residual or ambient light) and the irregular θ_d dependence of the spectrum S
  (a calibration difference between illumination angles is suspected) — the supplementary document of the paper may say.
- The shoulders of the lobe (θ_h 5°–20°) and the far field do not follow the factor model; they are dark, but a second,
  depolarized component would be the refinement.
- Edits keep S and N; physically, roughening a surface also depolarizes it.
- The source of the one-sided φ_d arms of the lobes: their direction differs between materials, so per-sample anisotropy or
  optics rather than one fixed instrument artifact. The symmetrization removes one-sided arms (aluminium) but only halves
  structure present at both φ_d and φ_d + 180° (SUJ2, Fake gold, black glass).
