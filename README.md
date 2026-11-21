# brdfkit

[![tests](https://github.com/Pero1031/brdfkit/actions/workflows/tests.yml/badge.svg)](https://github.com/Pero1031/brdfkit/actions/workflows/tests.yml)

Tools for measured BRDF tables (MERL, KAIST pBRDF 2020, hpBRDF 2025): loading, PCA compression for renderers,
and decomposition of the polarimetric data into physical quantities.
**Start with [docs/WALKTHROUGH.md](docs/WALKTHROUGH.md)**: from the downloaded hpBRDF tables to edited spectral, polarized materials, and how to repeat it.
Research code. The code is MIT; each dataset has its own terms ([LICENSES.md](LICENSES.md)).

## What it does

| Feature | Entry point | Description |
|---|---|---|
| Loading | `brdfkit.io.load(path)` | `.hpbrdf` / `.pbsdf` (Mitsuba tensor_file), `.mat` (KAIST 2020), `.binary` (MERL) into one `BRDFTable` |
| PCA compression | `brdfkit compress` / `eval` | A short latent vector per angular bin (hpBRDF: 13 GB → 231 MB, MERL: 10 MB). [docs/PCA.md](docs/PCA.md) |
| Export for renderers | `brdfkit export` | `.hpbpca` v2. One C++17 header, [`cpp/hpbrdf_pca.hpp`](cpp/hpbrdf_pca.hpp), reads and decodes it |
| Polarimetric analysis | `brdfkit lu-chipman` | Lu–Chipman decomposition (diattenuation, retardance, depolarization) compared with Fresnel and thin-film models. [docs/LU_CHIPMAN.md](docs/LU_CHIPMAN.md) |
| Inter-material PCA | `brdfkit merl-pca` | Matusik 2003-style PCA across the 100 MERL materials: common basis, material distribution, leave-one-out generalization error. [docs/MERL_PCA.md](docs/MERL_PCA.md) |
| Projection onto MERL PCA | `brdfkit merl-project` | Place a measured table (hpBRDF aluminium) in the MERL principal-component space: coefficients, residual, nearest materials. [docs/HPBRDF_ON_MERL.md](docs/HPBRDF_ON_MERL.md) |
| Instrument floor | `brdfkit floor`, `--floor` | Estimate the unpolarized offset that hpBRDF tables carry far from the lobe (on black glass) and remove it in `compress`, `eval`, `merl-project`. [docs/FLOOR.md](docs/FLOOR.md) |
| Calibration on black glass | `brdfkit calibrate`, `edit --calib` | Compare the black glass with Fresnel of a dielectric: the tables carry 0.47 (θ_d 20°) to 0.08 (75°) of the Fresnel energy, the same on aluminium; correct it, and see the spectral and polarimetric accuracy at the specular. [docs/CALIBRATION.md](docs/CALIBRATION.md) |
| Factorized model | `brdfkit factor` | m00 = grey angular shape A × spectrum S(θ_d, λ), polarization N(θ_d, λ): 6.8 MB, about 6% m00 error near the ridge; the parts can be recombined for editing. [docs/FACTOR.md](docs/FACTOR.md) |
| Latent editing | `brdfkit edit`, `examples/latent_demo.py` | Move the grey lobe of a factor model along the MERL PCA, keeping the measured spectrum and polarization: materials that were never measured. The lobe is first symmetrized over φ_d, which removes a streak the measured tables carry. [docs/LATENT_EDIT.md](docs/LATENT_EDIT.md) |
| Material swap | `brdfkit swap` | Keep the lobe shape of one material and take the colour (S), polarization (N) and reflectance of another, or of Fresnel from n, k: aluminium becomes gold, steel or black glass within 3% of the measured luminance and colour. [docs/SWAP.md](docs/SWAP.md) |
| Gloss knob | `brdfkit gloss` | Roughen a lobe to a target width (half width, or GGX-equivalent roughness) by convolving its microfacet slopes: one physical parameter that means the same for every material, colour and polarization kept. The measured lobes are the instrument's response (black glass, 1.1°) convolved with the surface. [docs/GLOSS.md](docs/GLOSS.md) |
| Quirks of a measured table | `brdfkit quirks` | One pass over a table, then the checks that matter before using it: inpainted angles, smoothing, the instrument's angular response, the floor, the radiometric scale against Fresnel, the spectrum along θ_d, one-sided φ_d arms, reciprocity, and the specular polarization against Fresnel. Results for eight hpBRDF materials in [docs/QUIRKS.md](docs/QUIRKS.md): what the tables hold is BRDF × cos(light) × g(θ_d), filed at φ_d + 180°, plus a floor |
| Which direction is the light | `brdfkit slots`, `--exchange` | The polarization of a diffuse body tells light from camera without a reference: the hpBRDF tables are filed at φ_d + 180° relative to the Mitsuba reader. [docs/SLOTS.md](docs/SLOTS.md) |
| Surface and body by polarization | `brdfkit umov` | On a coloured plastic the polarized part is the surface (white, Fresnel) and the rest the body: Umov's law over wavelength on one sample, and the surface reflection read out from under the body. [docs/UMOV.md](docs/UMOV.md) |
| Predicting an unmeasured plastic | `brdfkit predict` | A coloured plastic as a black one (surface) plus r(λ) × a white one's body, from one spectrum of the target; scored on two held-out plastics: the colour transfers (m00 within 4–6 % up to θ_d 55°), the surface's strength does not. [docs/PREDICT.md](docs/PREDICT.md) |
| The table as an ellipsometer | `brdfkit ellipso` | The phase Δ of the specular reflection is the part of a metal's polarization the tables get right; with the colour and the calibrated level it recognizes aluminium and steel among ten candidates, and leaves the fake gold unidentified (brass, TiN, gold and copper fit equally). [docs/ELLIPSO.md](docs/ELLIPSO.md) |
| Corrected tables for Mitsuba 3 | `brdfkit to-mitsuba` | A measured table or a factor model as a `tensor_file` for `measured_polarized`, with the floor removed, the directions exchanged, and the light's cosine and the glass's g(θ_d) divided out; checked against Mitsuba 3. [docs/MITSUBA.md](docs/MITSUBA.md) |
| GoLD_rayt integration | [`integrations/gold_rayt/`](integrations/gold_rayt/) | `Material` implementation and integration steps |

## Usage

```bash
pip install -e .            # numpy, scipy (add .[plot] for figures)
python tests/run_tests.py   # self-checks on synthetic data (and the real loaders when the files are present)

python -m brdfkit info     data/kaist/hpbrdf/Aluminum.hpbrdf
python -m brdfkit compress data/kaist/hpbrdf/Aluminum.hpbrdf out/pca/Aluminum_pca.npz --k1 8 --k2 32
python -m brdfkit eval     data/kaist/hpbrdf/Aluminum.hpbrdf out/pca/Aluminum_pca.npz --sweep
python -m brdfkit export   out/pca/Aluminum_pca.npz out/pca/Aluminum.hpbpca
python -m brdfkit lu-chipman out/pca/Aluminum_pca.npz --refidx data/kaist/refidx/Al_Rakic.yml --out out/analysis
python -m brdfkit merl-pca   data/merl --out out/merl_analysis/log
python -m brdfkit merl-project data/kaist/hpbrdf/Aluminum.hpbrdf --model out/merl_analysis/log/merl_pca.npz --out out/hpbrdf_on_merl --reference data/merl/aluminium.binary
python -m brdfkit to-mitsuba data/kaist/hpbrdf/Aluminum.hpbrdf out/clean/Aluminum.pbsdf --floor out/floor/Black_glass_floor.npz --exchange --cosine light --calib out/clean/black_glass_xc_cal.npz
```

Datasets are not included (see [LICENSES.md](LICENSES.md)). The examples and the optional real-data tests expect
`data/kaist/hpbrdf/Aluminum.hpbrdf` and the MERL `.binary` files in `data/merl/`; a symlink or junction into wherever
the datasets live works, or set `BRDFKIT_HPBRDF` / `BRDFKIT_MERL_DIR` for the tests.

From Python:

```python
from brdfkit.io import load
from brdfkit.decode import PCADecoder

t = load('data/merl/aluminium.binary')          # BRDFTable: axes, channels, values (or an iterator over phi_d slices)
dec = PCADecoder('out/pca/Aluminum_pca.npz')
r = dec.eval(to_light, to_view, 550.0)          # m00 (f_r without the cosine term), Mueller matrix M, reference vectors x_in / x_out
```

## Layout

```
brdfkit/         package: table (common representation), io/ (loaders), geometry (Rusinkiewicz coordinates, Stokes bases),
                 compress/pca, decode (reference decoder), export/hpbpca, analysis/{lu_chipman, merl_pca, merl_project, merl_subsets, floor, factor, latent_edit, calibrate, swap, gloss, quirks, slots, umov, clean, predict, ellipso},
                 color, cli
cpp/             home of the C++ decoder (hpbrdf_pca.hpp) and a driver that checks it against Python
integrations/    gold_rayt/: Material for GoLD_rayt, integration steps, build check
tests/           self-checks on synthetic BRDFs (no pytest needed: python tests/run_tests.py)
docs/            WALKTHROUGH.md (start here), QUIRKS.md, SLOTS.md, UMOV.md, PREDICT.md, ELLIPSO.md, MITSUBA.md, PCA.md, LU_CHIPMAN.md, MERL_PCA.md, HPBRDF_ON_MERL.md, FLOOR.md, FACTOR.md, LATENT_EDIT.md, CALIBRATION.md, SWAP.md, GLOSS.md
examples/        latent_demo.py: the whole chain from the downloaded data to rendered edits
data/, out/      datasets and derived files; not tracked by git
```

## Conventions

Directions are in the local shading frame (normal = +z): `to_light` (surface → light) and `to_view` (surface → viewer).
The Rusinkiewicz coordinates and the Mueller-matrix frames are those of the Mitsuba 3 `measured_polarized` plugin,
with reference vectors in the plane of reflection. `geometry.mueller_to_sp` rotates into the (s, p) bases of GoLD_rayt.
Details at the top of `brdfkit/geometry.py` and `cpp/hpbrdf_pca.hpp`.
