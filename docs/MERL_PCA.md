# MERL 100 材質の材質間 PCA（`brdfkit merl-pca`）

Matusik ら（2003）の「データ駆動反射モデル」の出発点になった解析を、MERL BRDF データベースの 100 材質で再現します。
各材質の BRDF 表全体を 1 本の特徴ベクトルにして材質をサンプルとする PCA を解き、共通基底と、
基底の中での材質の分布、材質を外して当てはめたときの汎化誤差（leave-one-out）を出します。
実装は `brdfkit/analysis/merl_pca.py`。

```bash
python -m brdfkit merl-pca D:\Jack_spallow\MERL_BRDF_database\brdfs --out out/merl_analysis/log
python -m brdfkit merl-pca D:\Jack_spallow\MERL_BRDF_database\brdfs --out out/merl_analysis/nielsen --mapping nielsen
```

| オプション | 意味 | 既定 |
|---|---|---|
| `--stride` | 格子の間引き。1 つの整数か `phi_d,theta_d,theta_h` | 2（90 × 45 × 45 ビン、材質あたり 1.7 MB） |
| `--k` | 保存する基底の本数（係数と寄与率は全成分を保存） | 16 |
| `--mapping` | 特徴量の変換。`log` = log(max(f, ε))（Matusik 2003）、`nielsen` = log((f cos θi cos θo + ε) / (材質中央値 + ε))（Nielsen, Jensen, Ramamoorthi 2015） | log |
| `--eps` | ε | log: 1e-5、nielsen: 1e-3 |
| `--no-plot` | 図を出さない | |

## 方法

- 全材質で測定済みのビン（MERL では幾何的に決まる 79% が全材質共通）× RGB 3 チャンネルを並べて 1 材質 = 1 ベクトル（stride 2 で D ≈ 44 万次元）。
- D × D の共分散は作らず、中心化したベクトルの Gram 行列（100 × 100）を固有分解する。全固有値（= 全成分の寄与率）、
  全材質の全成分の係数、各材質の「k 本で切ったときの残差」がこの 100 × 100 行列だけから厳密に出る。
- leave-one-out: 材質 i を外した 99 材質で PCA を作り、i をその基底に射影したときの残差。これも Gram 行列の代数だけで出る
  （99 × 99 の固有分解を 100 回）。再構成誤差との差が、基底が「その材質を見ていたから」合っている分。
- 基底は先頭 k 本だけ D 次元に戻して保存する。`MaterialPCA.project(values)` で同じ格子上の別の表を係数にでき、
  `reconstruct(coef)` で係数から（間引いた格子上の）線形 BRDF に戻せる。

## 出力（`--out` 以下）

- `merl_pca.npz`: `names`、`mean [D]`、`basis [D, k]`、`coef [n, n]`（全成分）、`eigvals`、`evr`、`valid [P', T', H']`、間引いた角度軸、
  `mapping`、`eps`、`ref`（nielsen の中央値）、`cosw`（cos θi cos θo）、`recon_rmse [n, n+1]`、`loo_rmse [n, n−1]`、`meta`。
  `MaterialPCA.load(path)` で読む。
- `merl_pca_report.txt`: k ごとの累積寄与率と誤差、材質ごとの係数・誤差・最近傍材質。
- `merl_pca_variance.png`（寄与率と誤差の k 依存）、`merl_pca_scatter.png`（PC1–PC2、PC3–PC4 の材質散布図。色は鏡面ピーク / 中央値）、
  `merl_pca_basis.png`（平均と PC1〜6 の θh × θd 像、φd = 90°）、`merl_pca_recon.png`（LOO 誤差が最小の 1 材質と最大の 3 材質の再構成）。

## 結果（100 材質、stride 2、2026-11-16）

| k | 累積寄与率 log / nielsen | log-RMSE 再構成 中央値 | LOO 中央値 / p90 / 最大 | nielsen: 再構成 / LOO 中央値 |
|---|---|---|---|---|
| 1 | 0.588 / 0.611 | 1.01 | 1.03 / 1.51 / 2.68 | 0.84 / 0.86 |
| 2 | 0.814 / 0.814 | 0.64 | 0.66 / 1.15 / 1.35 | 0.55 / 0.57 |
| 3 | 0.893 / 0.894 | 0.49 | 0.51 / 0.82 / 1.19 | 0.44 / 0.46 |
| 5 | 0.965 / 0.964 | 0.30 | 0.32 / 0.47 / 0.67 | 0.25 / 0.28 |
| 8 | 0.980 / 0.981 | 0.22 | 0.24 / 0.38 / 0.57 | 0.18 / 0.20 |
| 16 | 0.993 / 0.993 | 0.13 | 0.17 / 0.28 / 0.46 | 0.11 / 0.13 |
| 32 | 0.997 / 0.998 | 0.09 | 0.13 / 0.25 / 0.42 | 0.07 / 0.10 |
| 64 | 0.9995 / 0.9996 | 0.04 | 0.10 / 0.21 / 0.39 | 0.03 / 0.08 |

log-RMSE は自然対数の単位で、0.1 が「ビンあたり約 10% の誤差」。所要時間は読み込み 41 s、PCA（LOO 込み）97 s。

**寄与率と誤差の乖離。** 5 成分で分散の 96% を説明するが、ビンあたりの対数誤差はまだ 30%（中央値）ある。
分散は「金属か拡散か」（PC1、59%）と「暗いか明るいか」（PC2、22%）の大きな対比に支配されていて、
鏡面ローブの細かい形は残り数% の中に入っているため。レンダリングで使える精度（対数で 0.1 程度）には 16〜32 成分が要る。

**汎化。** 材質を外して当てはめると、k = 16 で中央値 0.17、最悪 0.46（grease-covered-steel）。
再構成との差は k = 16 で中央値 0.03 だが、金属では 0.1〜0.3 開く（steel 0.22 → 0.40、chrome-steel 0.27 → 0.43）。
100 材質の中で金属は十数個しかなく、鏡面ピークの形の多様性を基底が十分に持っていない。
LOO が悪い順: grease-covered-steel、chrome-steel、steel、black-obsidian、green-acrylic、chrome、aluminium、ss440、tungsten-carbide。
良い順: pink-fabric（0.07）、gold-paint、white-fabric、neoprene-rubber、pink-felt、polyurethane-foam。布や塗料は基底でよく表せる。

**主成分の意味**（`merl_pca_basis.png`、`merl_pca_scatter.png`）。

- PC1（59%）: 鏡面 vs 拡散。θh < 10° を下げ、それ以外を上げる。金属（chrome −1900、tungsten-carbide、steel、aluminium −1170）が負、
  布（pink-fabric、white-fabric、polyethylene +1190）が正。
- PC2（22%）: 全体の明るさと掠め角。black-obsidian、black-fabric、black-phenolic が正、alumina-oxide、white-acrylic、white-marble が負。
- PC3（8%）: ローブの幅。θh 10〜30° を下げて 50° 以上を上げる。specular-*-phenolic、maroon-plastic が正、metallic paint と nickel が負。
- PC4（5%）: 色相。blue-acrylic、green-plastic、specular-blue-phenolic が正、red-metallic-paint、specular-red/orange-phenolic が負（図はチャンネル平均なので像には出ない）。
- PC5 以降は 2% 以下。PC6 で 0.6%。

**最近傍。** 16 次元の係数空間での最近傍は、名前の対（red-fabric ↔ red-fabric2、silver-paint ↔ pearl-paint、cherry-235 ↔ ipswich-pine-221、
gold-metallic-paint2 ↔ two-layer-gold、chrome ↔ tungsten-carbide）がそのまま出る。aluminium の最近傍は brass、black-obsidian は specular-black-phenolic。

**nielsen 写像。** cos θi cos θo を掛けて材質の中央値で割ると、全 k で誤差が 15〜25% 小さい（k = 16 の再構成 0.11、LOO 0.13）。
掠め角の大きな値が cos で抑えられ、材質共通の形が中央値で消えて、PCA が材質間の差だけを説明すればよくなるため。
一方で cos θi cos θo = 0 の境界ビンが使えなくなる（共通ビン 77%）。

**線形空間の誤差。** 当てはめた材質を k 本で再構成したときの f cos θi cos θo の相対 RMSE は k = 16 で中央値 0.17、p90 0.32。
鏡面ピークの直近が支配するので、`docs/PCA.md` の角度ビンごとの PCA（m00 相対 RMSE 0.2%）とは目的も精度も違う。
こちらは「材質の空間」を見るための道具で、1 材質をレンダラーに渡す圧縮には角度ビンごとの PCA を使う。

## hpBRDF との接続

`MaterialPCA.project` は MERL と同じ格子（180 × 90 × 90、RGB）の表を受け取る。hpBRDF の Aluminum（68 波長、ミュラー行列、91 × 91 × 361 の格子）を
この空間に置くには、m00 を RGB に落として MERL 格子に再サンプルする変換が要る。これは未実装。
