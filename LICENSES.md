# Dataset licenses

The code in this repository is MIT (`LICENSE`); the datasets it reads have their own, separate terms.
**Neither the datasets nor the compressed tables derived from them (`.npz`, `.hpbpca`) are committed** (`data/` and `out/` are in `.gitignore`).
Cite the respective papers when using them in publications.

| Dataset | Source | Terms | Citation |
|---|---|---|---|
| hpBRDF (2025) | [Hugging Face](https://huggingface.co/datasets/yunseongmoon/Hyperspectral-Polarimetric-BRDF), [GitHub](https://github.com/yunseong0518/hpBRDF) | **No license stated** (no dataset card, no LICENSE file). Reading it for research is what the release is for; ask the authors (Seung-Hwan Baek, POSTECH) before redistributing or using it commercially | Moon, Maeda, Shin, Hwang, Kim, Kim, Baek, "Hyperspectral Polarimetric BRDFs of Real-world Materials", SIGGRAPH Asia 2025 |
| KAIST pBRDF (2020) | [Distribution page](https://vclab.kaist.ac.kr/siggraph2020/pbrdfdataset/kaistdataset.html) | No formal license; only "cite the paper if you use it" | Baek, Zeltner, Ku, Hwang, Tong, Jakob, Kim, "Image-Based Acquisition and Modeling of Polarimetric Reflectance", ACM TOG 39(4), SIGGRAPH 2020 |
| MERL BRDF (2006 package) | merl.com (historical) | Use, copying and modification allowed for educational, research and non-profit purposes; keep the three-paragraph copyright notice and disclaimer on every copy. Commercial use needs permission from MERL (license@merl.com). The bundled `BRDFRead.cpp` is MERL's work (brdfkit's reader is an independent implementation) | Matusik, Pfister, Brand, McMillan, "A Data-Driven Reflectance Model", ACM TOG 22(3), 2003 |
| MERL BRDF (2023 Zenodo redistribution) | [doi:10.5281/zenodo.8101680](https://doi.org/10.5281/zenodo.8101680) | **CC BY-SA 4.0** (attribution, share-alike) | same as above |
| Refractive-index data (`Al_Rakic.yml` etc.) | [refractiveindex.info](https://refractiveindex.info) | Public domain (CC0) | reference at the top of each file |

Read as reference implementations (no code copied):
the Mitsuba 3 `measured_polarized` plugin (BSD-3) and `pbrdf_io.py` from the hpBRDF repository.
