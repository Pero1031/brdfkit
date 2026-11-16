# brdfkit

[![tests](https://github.com/Pero1031/brdfkit/actions/workflows/tests.yml/badge.svg)](https://github.com/Pero1031/brdfkit/actions/workflows/tests.yml)

Tools for measured BRDF tables (MERL, KAIST pBRDF 2020, hpBRDF 2025): loading, PCA compression for renderers,
and decomposition of the polarimetric data into physical quantities.
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
                 compress/pca, decode (reference decoder), export/hpbpca, analysis/{lu_chipman, merl_pca, merl_project}, color, cli
cpp/             home of the C++ decoder (hpbrdf_pca.hpp) and a driver that checks it against Python
integrations/    gold_rayt/: Material for GoLD_rayt, integration steps, build check
tests/           self-checks on synthetic BRDFs (no pytest needed: python tests/run_tests.py)
docs/            PCA.md, LU_CHIPMAN.md, MERL_PCA.md, HPBRDF_ON_MERL.md
data/, out/      datasets and derived files; not tracked by git
```

## Conventions

Directions are in the local shading frame (normal = +z): `to_light` (surface → light) and `to_view` (surface → viewer).
The Rusinkiewicz coordinates and the Mueller-matrix frames are those of the Mitsuba 3 `measured_polarized` plugin,
with reference vectors in the plane of reflection. `geometry.mueller_to_sp` rotates into the (s, p) bases of GoLD_rayt.
Details at the top of `brdfkit/geometry.py` and `cpp/hpbrdf_pca.hpp`.
