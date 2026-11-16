# GoLD_rayt への組み込み手順（コード差分）

ヘッダ 2 枚を置き、シーンファイルの読み込みに `"type": "measured"` を 3 か所で足すだけです。新しい .cpp はありません。

## 1. ファイルを置く

| バンドル内 | GoLD_rayt 内の置き場所 |
|---|---|
| `include/Materials/hpbrdf_pca.hpp` | `GoLD_rayt/include/Materials/hpbrdf_pca.hpp` |
| `include/Materials/MeasuredHpBRDF.hpp` | `GoLD_rayt/include/Materials/MeasuredHpBRDF.hpp` |
| `out/pca/Aluminum.hpbpca`（brdfkit の出力） | `GoLD_rayt/assets/hpbrdf/Aluminum.hpbpca`（231 MB、.gitignore に入れる） |
| `scenes/hpbrdf_aluminum.json` | `GoLD_rayt/scenes/hpbrdf_aluminum.json` |

`GoLD_rayt.vcxproj` の `<ClInclude>` に 2 つのヘッダを足すとソリューションエクスプローラーに出ます（ビルドには不要）。
`.filters` では `Materials` フィルタに入れてください。

## 2. `include/Cli/SceneDescription.hpp` — MaterialSpec

```cpp
struct MaterialSpec {
    enum class Type { Diffuse, Conductor, Dielectric, Emission, ThinFilm, Layered, Hair, Skin, Measured };   // Measured を追加
    ...
    // measured: a PCA-compressed hpBRDF table (Materials/MeasuredHpBRDF.hpp)
    std::string measuredFile;                      //   the .hpbpca file
    double alphaSample = 0.1;                      //   GGX alpha of the sampling lobe (only the noise depends on it)
```

## 3. `src/Cli/SceneFileMaterials.cpp` — ParseMaterial

`type == "emission"` の分岐の手前あたりに追加:

```cpp
    else if (type == "measured") {   // { "type": "measured", "file": "assets/hpbrdf/Aluminum.hpbpca", "alphaSample": 0.1 }
        CheckKeys(m, { "type", "file", "alphaSample" }, where);
        s.type = MaterialSpec::Type::Measured;
        if (!m.contains("file")) Fail(where, "needs a \"file\": a PCA-compressed hpBRDF table (.hpbpca)");
        s.measuredFile = ReadString(m.at("file"), where + ".file");
        Optional<double>(m, "alphaSample", where, s.alphaSample, ReadNumber);
        if (!(s.alphaSample > 0)) Fail(where + ".alphaSample", "must be above 0");
    }
```

先頭のエラーメッセージ（`"must be an object with a \"type\" (diffuse, conductor, ..."`）と、
最後の `else Fail(...)` の型一覧にも `measured` を足してください。

## 4. `src/Cli/Scenes/FileScene.cpp` — 材質の生成

include を追加:

```cpp
#include "Materials/MeasuredHpBRDF.hpp"
```

`switch (m.type)` に case を追加（`case MaterialSpec::Type::Conductor` の隣）:

```cpp
        case MaterialSpec::Type::Measured: {
            std::string err;
            auto mat = MeasuredHpBRDF::load(m.measuredFile, Real(m.alphaSample), &err);
            if (!mat) {
                std::cerr << "[SceneFile] material \"" << name << "\": " << m.measuredFile << ": " << err << "\n";
                return 1;
            }
            const auto& t = mat->table();
            std::cout << "[SceneFile] material \"" << name << "\": measured hpBRDF " << m.measuredFile
                << " (" << t.n_phi << "x" << t.n_thd << "x" << t.n_thh << " bins, " << t.n_wvl << " bands "
                << t.wvls.front() << "-" << t.wvls.back() << " nm, latent " << t.k1 << "+" << t.k2 << ")\n";
            materials[name] = mat;
            break;
        }
```

`MaterialSpec::Type` で switch している場所は `FileScene.cpp` の 1 か所だけです（`grep -rn 'case MaterialSpec::Type'`）。
`/W4` で列挙子の漏れ警告が出たら、そこに `Measured` を足してください。

## 5. `docs/SCENE_FILE.md` の材質の表に 1 行

```
| `measured` | `file`（PCA 圧縮した hpBRDF の表 `.hpbpca`。`python -m brdfkit export` で作る）、`alphaSample`（重点サンプリングに使う GGX の α。絵は変わらずノイズだけ変わる。既定 0.1） |
```

## 6. 動作確認

```bat
GoLD_rayt.exe scenes\hpbrdf_aluminum.json --spp 16 --width 400
GoLD_rayt.exe scenes\hpbrdf_aluminum.json --stokes --polarizer 0,90
```

- 読み込み時に `[SceneFile] material "aluminum": measured hpBRDF ...` が出る。RGB ビルドは RGB 表の焼き込みに数秒かかります（OpenMP 並列）。
- `--stokes` の DoLP 画像で、アルミ球の縁（斜入射）ほど偏光度が高く、正面でほぼ 0 になれば、ミュラー行列の受け渡しは正しく動いています（金球と同じ傾向）。
- 偏光フィルター 0° と 90° で、球の上下（入射面が縦）と左右で明暗が入れ替わります。

## 7. スペクトルビルド

`Core/RenderMode.hpp` の `RAYT_ENABLE_SPECTRAL 1` でビルドし直すと、eval が経路の 4 波長それぞれで表を復号します。
GoLD_rayt の波長範囲 360〜830 nm のうち 360〜414 nm は表の外なので 414 nm の値で延ばします（アルミはほぼ平坦なので影響は小さい）。
