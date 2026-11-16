# brdfkit

実測 BRDF 表（MERL、KAIST pBRDF 2020、hpBRDF 2025）を読み込み、PCA で圧縮してレンダラーに渡し、偏光の物理量に分解するためのツール群です。
研究用。コードは MIT、データセットは各自の条件に従います（[LICENSES.md](LICENSES.md)）。

## できること

| 機能 | 入口 | 説明 |
|---|---|---|
| 読み込み | `brdfkit.io.load(path)` | `.hpbrdf` / `.pbsdf`（Mitsuba tensor_file）、`.mat`（KAIST 2020）、`.binary`（MERL）を共通の `BRDFTable` に |
| PCA 圧縮 | `brdfkit compress` / `eval` | 角度ビンごとの潜在ベクトル（hpBRDF: 13 GB → 231 MB、MERL: 10 MB）。[docs/PCA.md](docs/PCA.md) |
| レンダラー向け書き出し | `brdfkit export` | `.hpbpca` v2。C++17 のヘッダ 1 枚 [`cpp/hpbrdf_pca.hpp`](cpp/hpbrdf_pca.hpp) で読み込みから復号まで |
| 偏光分析 | `brdfkit lu-chipman` | Lu–Chipman 分解（二色性・位相差・偏光解消）と Fresnel / 薄膜モデルとの比較。[docs/LU_CHIPMAN.md](docs/LU_CHIPMAN.md) |
| GoLD_rayt 連携 | [`integrations/gold_rayt/`](integrations/gold_rayt/) | `Material` 実装と組み込み手順 |

## 使い方

```bash
pip install -e .            # numpy, scipy（図を出すなら .[plot]）
python tests/run_tests.py   # 合成データでの自己検証（実データがあればその読み込みも）

python -m brdfkit info     data/kaist/hpbrdf/Aluminum.hpbrdf
python -m brdfkit compress data/kaist/hpbrdf/Aluminum.hpbrdf out/pca/Aluminum_pca.npz --k1 8 --k2 32
python -m brdfkit eval     data/kaist/hpbrdf/Aluminum.hpbrdf out/pca/Aluminum_pca.npz --sweep
python -m brdfkit export   out/pca/Aluminum_pca.npz out/pca/Aluminum.hpbpca
python -m brdfkit lu-chipman out/pca/Aluminum_pca.npz --refidx data/kaist/refidx/Al_Rakic.yml --out out/analysis
```

Python から:

```python
from brdfkit.io import load
from brdfkit.decode import PCADecoder

t = load('data/merl/aluminium.binary')          # BRDFTable: 軸、チャンネル、値（または φ_d スライスのイテレータ）
dec = PCADecoder('out/pca/Aluminum_pca.npz')
r = dec.eval(to_light, to_view, 550.0)          # m00（f_r、コサイン項なし）とミュラー行列 M、基準ベクトル x_in / x_out
```

## 構成

```
brdfkit/         パッケージ: table（共通表現）、io/（loaders）、geometry（Rusinkiewicz 座標、Stokes 基底）、
                 compress/pca、decode（参照復号器）、export/hpbpca、analysis/lu_chipman、cli
cpp/             C++ 復号器の本籍（hpbrdf_pca.hpp）と Python との一致確認ドライバ
integrations/    gold_rayt/: GoLD_rayt 用 Material、組み込み手順、ビルド確認
tests/           合成 BRDF での自己検証（pytest 不要: python tests/run_tests.py）
docs/            PCA.md、LU_CHIPMAN.md
data/, out/      データセットと派生物。git 対象外
```

## 規約

方向はローカルシェーディング座標（法線 +z）で、`to_light`（表面→光源）と `to_view`（表面→視点）。
Rusinkiewicz 座標とミュラー行列のフレームは Mitsuba 3 の `measured_polarized` プラグインと同じで、
基準ベクトルは反射面内にあります。GoLD_rayt の (s, p) 基底へは `geometry.mueller_to_sp` で回します。
詳細は `brdfkit/geometry.py` と `cpp/hpbrdf_pca.hpp` の冒頭。
