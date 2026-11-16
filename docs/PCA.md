# PCA 圧縮（`brdfkit compress` / `eval` / `export`）

実測 BRDF 表を、角度ビンごとの短い潜在ベクトルと小さな基底に圧縮します。
13 GB の hpBRDF をメモリに載せず、φ_d スライス（約 36 MB）単位でストリーミング処理します。

```bash
python -m brdfkit info     data/kaist/hpbrdf/Aluminum.hpbrdf
python -m brdfkit compress data/kaist/hpbrdf/Aluminum.hpbrdf out/pca/Aluminum_pca.npz --k1 8 --k2 32
python -m brdfkit eval     data/kaist/hpbrdf/Aluminum.hpbrdf out/pca/Aluminum_pca.npz --sweep
python -m brdfkit export   out/pca/Aluminum_pca.npz out/pca/Aluminum.hpbpca
```

| オプション | 意味 | 既定 |
|---|---|---|
| `--k1` | log m00 のチャンネル方向（68 波長、または RGB 3）の主成分数 | 8 |
| `--k2` | 正規化ミュラー成分 m_ij/m00（チャンネル × 15）の主成分数。無偏光の表では 0 | 32 |
| `--dtype` | 係数の保存型（`float16` / `float32`） | float16 |
| `--step` | 共分散推定時の φ_d ストライド（高速化用） | 1 |

## 表現

1 つの角度ビンを 2 つの PCA に分けて表します。

- **Part A（強度）** `x1 = log(max(m00, eps))`、チャンネル数次元 → k1 係数。
  対数にすることで鏡面ピークのダイナミックレンジを抑え、PCA がピークだけを説明するのを防ぎます。
- **Part B（偏光、偏光表のみ）** `x2 = clip(m_ij / m00, −1, 1)`（ij ≠ 00）、チャンネル × 15 次元 → k2 係数。
  受動的なミュラー行列は |m_ij| ≤ m00 なので、[−1, 1] へのクリップは物理的に正当です。

復号:

```
m00(ch)   = exp( mean_logm00 + basis_logm00 · c1 )
m_ij(ch)  = m00(ch) · ( mean_ratio + basis_ratio · c2 )[ratio_index]
```

未測定ビン（NaN）は学習と評価から外し、係数表では同じ φ_d スライス内の最近傍有効ビンの係数で埋めます。
`valid` マスクは表に残すので、レンダラー側で別の埋め方に差し替えられます。

## 出力

- `*.npz`: 係数 `coef_logm00` `coef_ratio`、`valid`、平均と基底、角度軸、`channels`、`channel_kind`、`polarized`、`meta`（JSON）。
  復号の参照実装は `brdfkit/decode.py`（`PCADecoder`）。
- `*_basis_full.npz`: 全固有ベクトルと寄与率（`eval --sweep` 用）。
- `*.hpbpca`（`export`）: レンダラー向けの平坦なバイナリ（形式 v2、レイアウトは `cpp/hpbrdf_pca.hpp` 冒頭と `brdfkit/export/hpbpca.py`）。
  係数は 1 ビンあたり k1 + k2 個を連続で置き、`cpp/hpbrdf_pca.hpp` が読み込みから復号まで行います。

## 結果

**hpBRDF Aluminum（k1 = 8、k2 = 32、float16）**

| 項目 | 値 |
|---|---|
| 元ファイル | 13.0 GB |
| 圧縮後 `.npz` / `.hpbpca` | 231 MB |
| 処理時間 | 共分散 209 s + 符号化 105 s |
| 累積寄与率 | log m00: 8 成分で 99.85%、ratio: 32 成分で 99.4% |
| m00 各波長の相対 RMSE | 10.3%（ほぼ測定ノイズ。残差は波長方向に白色で、2 階差分から見積もったノイズ水準と一致） |
| m00 を 68 バンドで積分した値の相対誤差 | 中央値 0.2〜0.4%、99 パーセンタイル 約 5% |
| 正規化ミュラー成分の RMSE | 0.007〜0.023 |
| 法線入射の方向アルベド | 元テーブル 0.353、復号後 0.351 |

**MERL aluminium（RGB、k1 = 3、float16）**: 9.7 MB、m00 相対 RMSE 0.23%（3 チャンネルを 3 成分で持つので float16 の丸めだけ）。

**既知の制限**: 鏡面ピークの直近では元テーブルより 10〜18% 低くなることがあります（係数空間の補間 ＝ log 空間の補間）。
ピークの精度が要るなら θ_h の解像度を上げるか、ピーク周辺だけ元テーブルを参照する構成になります。
