# データセットのライセンス

このリポジトリのコードは MIT（`LICENSE`）ですが、扱うデータセットには別々の条件があります。
**データ本体と、そこから作った圧縮テーブル（`.npz`、`.hpbpca`）はコミットしません**（`data/`、`out/` は `.gitignore`）。
論文や公開物で使うときは、それぞれの引用を付けてください。

| データセット | 入手先 | 条件 | 引用 |
|---|---|---|---|
| hpBRDF（2025） | [Hugging Face](https://huggingface.co/datasets/yunseongmoon/Hyperspectral-Polarimetric-BRDF)、[GitHub](https://github.com/yunseong0518/hpBRDF) | **ライセンス表記なし**（dataset card なし、LICENSE なし）。研究目的での読み込みは配布の意図通りだが、再配布・商用は著者（Seung-Hwan Baek, POSTECH）に確認する | Moon, Maeda, Shin, Hwang, Kim, Kim, Baek, "Hyperspectral Polarimetric BRDFs of Real-world Materials", SIGGRAPH Asia 2025 |
| KAIST pBRDF（2020） | [配布ページ](https://vclab.kaist.ac.kr/siggraph2020/pbrdfdataset/kaistdataset.html) | 正式なライセンスなし。「使ったら論文を引用」とだけ | Baek, Zeltner, Ku, Hwang, Tong, Jakob, Kim, "Image-Based Acquisition and Modeling of Polarimetric Reflectance", ACM TOG 39(4), SIGGRAPH 2020 |
| MERL BRDF（2006 パッケージ） | merl.com（旧） | 教育・研究・非営利目的なら使用・複製・改変可。著作権表示と免責の 3 段落を全コピーに残す。商用は MERL（license@merl.com）に許諾を求める。同梱の `BRDFRead.cpp` は MERL の著作物（brdfkit の読み込みは独自実装） | Matusik, Pfister, Brand, McMillan, "A Data-Driven Reflectance Model", ACM TOG 22(3), 2003 |
| MERL BRDF（2023 Zenodo 再配布） | [doi:10.5281/zenodo.8101680](https://doi.org/10.5281/zenodo.8101680) | **CC BY-SA 4.0**（表示・継承） | 同上 |
| 屈折率データ（`Al_Rakic.yml` など） | [refractiveindex.info](https://refractiveindex.info) | パブリックドメイン（CC0） | 各ファイル冒頭の文献 |

参照実装として読んだもの（コードはコピーしていません）:
Mitsuba 3 の `measured_polarized` プラグイン（BSD-3）、hpBRDF リポジトリの `pbrdf_io.py`。
