# hpBRDF（と brdfkit が書き出す他の BRDF 表）を GoLD_rayt で描くための移植バンドル

KAIST/POSTECH の hpBRDF（実測の分光・偏光 BRDF、1 材質 13 GB）を PCA で 231 MB に圧縮し、
GoLD_rayt の `Material` として評価・サンプリング・偏光（`scatterMueller`）まで動く形にしたものです。
GoLD_rayt 側のコードを読んで API と規約に合わせて書き、GoLD_rayt のヘッダと `GGX.cpp` に対してビルド・リンク・実行まで確認しました
（`test/build_check.bat`、RGB ビルドとスペクトルビルドの両方）。GoLD_rayt 本体に組み込んでの描画はまだです。
手順は [PATCH_NOTES.md](PATCH_NOTES.md) にあります。

## 中身

| ファイル | 役割 |
|---|---|
| `include/Materials/MeasuredHpBRDF.hpp` | GoLD_rayt の `Material` 実装（ヘッダのみ）。eval / sample / pdf / scatterMueller / aovAlbedo |
| `include/Materials/hpbrdf_pca.hpp` | 依存なしの復号器（`namespace hpbrdf`）。ファイル読み込み、Rusinkiewicz 座標、潜在ベクトルの補間、復号 |
| `out/pca/Aluminum.hpbpca`（リポジトリ直下の out/、git 対象外） | 圧縮済みアルミ（231 MB、潜在 8+32 次元、float16）。`python -m brdfkit export` で作る |
| `scenes/hpbrdf_aluminum.json` | 例: 実測アルミ球と Johnson & Christy の金球を並べたシーン |
| `test/` | GoLD_rayt のヘッダに対してコンパイルが通ることだけを確かめた TU とバッチ |

## GoLD_rayt の規約への対応

| GoLD_rayt | 復号器（Mitsuba 3 の `measured_polarized` 規約） |
|---|---|
| `wo`（視点へ）、`wi`（光源へ）、ワールド座標 | `to_view`、`to_light`、`frame::Frame(rec.n)` のローカル座標に変換して渡す |
| `eval` はコサイン項なしの f_r | 表の m00 もコサイン項なし。そのまま返す |
| Spectrum = 4 波長（スペクトルビルド） | 波長ごとに `m00(c, λ)` を復号 |
| Spectrum = RGB（RGB ビルド） | 読み込み時に 68 バンドの m00 を等色関数で RGB に畳んだ表を焼く（約 3 M ビン × 3）。等エネルギー白で正規化 |
| `scatterMueller`: M/M00 を s 基底（s ⟂ 入射面、両ビーム共通）で | 表のフレーム（基準ベクトルが反射面内）から `pol::basisChange` で回す: `basisChange(x_out, s, wo) · M/m00 · basisChange(s, x_in, −wi)` |
| Stokes: (x, y = k × x, k)、exp(−iωt) | Mitsuba 3 も右手系 (x, y, k)。符号は下記で照合済み |
| 重点サンプリング | GGX VNDF（`alphaSample`）0.9 + コサイン 0.1 の混合。どの α でも結果は不偏 |

## 偏光の符号規約の照合

テーブルのミュラー行列を上の式で GoLD の s 基底に回し、鏡面反射の配置（θ_h = 0）で、
GoLD 自身の `pol::fresnelAmplitudes` + `setInterface` が Al（Rakić の n, k）に対して出すミュラー行列と比べました。

| 要素 | 表（450〜800 nm、θ = 20〜75°） | GoLD の Fresnel | 符号 |
|---|---|---|---|
| m01 | +0.06 〜 +0.48 | +0.005 〜 +0.22 | 一致（s 偏光が強く反射） |
| m22, m33 | 負（掠め角で 0 を越える） | 負 | 一致 |
| m23 | +0.06 〜 +0.91 | +0.03 〜 +0.91 | 一致 |
| m32 | −0.02 〜 −0.79 | −0.03 〜 −0.91 | 一致 |

全要素で符号が一致したので、S3（円偏光）の符号反転は入れていません。
大きさの差は、実測のアルミ表面が理想の Fresnel 面ではない（酸化膜・粗さで m01 が大きく、|m22| が小さい）ためで、
これこそが実測 BRDF を使う意味です。

## 動作確認（`test/check_rgb.exe`、`test/check_spectral.exe`、バンドル直下で実行）

| 項目 | RGB ビルド | スペクトルビルド（450/550/650/800 nm） |
|---|---|---|
| 読み込み（RGB 表の焼き込み込み） | 0.85 秒 | 0.35 秒 |
| f_r 鏡面方向（θ = 45°） | 11.8 / 13.7 / 12.7 | 12.2 / 13.7 / 11.6 / 11.3 |
| f_r 鏡面から外れた方向 | 0.018 / 0.021 / 0.034 | 0.033 / 0.021 / 0.019 / 0.016 |
| `scatterMueller`（鏡面、450 nm）m01 / m22 / m23 / m32 / m33 | — | +0.127 / −0.792 / +0.428 / −0.396 / −0.874（Python 側の照合値と一致） |
| サンプリング E[f cos / pdf]（θ_o = 45°、2 万本） | 0.231 | 0.267（450 nm） |
| eval 1 回 | 441 ns | 776 ns |

**アルミの明るさについて。** この実測アルミの方向アルベドは法線入射で 0.35、45° で 0.26、75° で 0.08 と、
理想の Al の Fresnel 反射率（約 0.9）よりずっと低く、斜めほど暗くなります。圧縮で失われたのではなく
（元テーブルを直接積分しても 0.353、復号後 0.351）、データセットの試料がそういう反射率です。
`conductor` の Au と並べると暗く見えますが、それが実測値です。

**相反性。** 表は厳密には相反ではありません（wo と wi を入れ替えると鏡面外で 1〜2 割違う所があります）。
Mitsuba のプラグインと同じく、そのまま使っています。必要なら圧縮時に φ_d と φ_d + π を平均して対称化できます。

## 性能とメモリ

| 項目 | 値 |
|---|---|
| 読み込み | 0.33 秒（float16 → float32 展開）＋ RGB ビルドは RGB 表の焼き込み数秒 |
| RAM | 係数 456 MB（float32）＋ RGB 表 36 MB |
| 1 回の評価 | 無偏光 約 0.65 µs、偏光 約 1.0 µs（単体ベンチ、GoLD_rayt の中では未計測） |

メモリが気になるなら、係数を float16 のまま持って補間時に変換する形に `hpbrdf_pca.hpp` を変えられます（約 230 MB に半減）。

## 既知の制限

- 掠め角の約 25% のビンは未測定で、同じ φ_d スライス内の最近傍ビンの値が入っています。`valid` マスクは表に残してあります。
- 鏡面ピークの直近では元テーブルより 10〜18% 低くなることがあります（係数空間の補間 ＝ log 空間の補間）。
- 波長 414〜950 nm の外は端のバンドで延ばします。GoLD_rayt の 360〜414 nm がこれに当たります。
- RGB ビルドの `scatterMueller` は R/G/B を 610/550/465 nm で代表させています。
- 他の 13 材質（黒ガラス、各色の粗いプラスチック、シリコーン、SUJ2、Fake gold、ビリヤード球）は
  Hugging Face から取得して `python -m brdfkit compress` → `export` で同じ形式にできます（1 材質 約 5 分）。
