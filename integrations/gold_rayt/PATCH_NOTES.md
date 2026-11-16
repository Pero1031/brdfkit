# Integration into GoLD_rayt (code changes)

Drop in two headers and add `"type": "measured"` to the scene-file loader in three places. No new .cpp files.

## 1. Place the files

| In this bundle | Where it goes in GoLD_rayt |
|---|---|
| `include/Materials/hpbrdf_pca.hpp` | `GoLD_rayt/include/Materials/hpbrdf_pca.hpp` |
| `include/Materials/MeasuredHpBRDF.hpp` | `GoLD_rayt/include/Materials/MeasuredHpBRDF.hpp` |
| `out/pca/Aluminum.hpbpca` (brdfkit output) | `GoLD_rayt/assets/hpbrdf/Aluminum.hpbpca` (231 MB; add to .gitignore) |
| `scenes/hpbrdf_aluminum.json` | `GoLD_rayt/scenes/hpbrdf_aluminum.json` |

Adding the two headers to `<ClInclude>` in `GoLD_rayt.vcxproj` shows them in the Solution Explorer (not needed for the build).
In `.filters`, put them under the `Materials` filter.

## 2. `include/Cli/SceneDescription.hpp` — MaterialSpec

```cpp
struct MaterialSpec {
    enum class Type { Diffuse, Conductor, Dielectric, Emission, ThinFilm, Layered, Hair, Skin, Measured };   // add Measured
    ...
    // measured: a PCA-compressed hpBRDF table (Materials/MeasuredHpBRDF.hpp)
    std::string measuredFile;                      //   the .hpbpca file
    double alphaSample = 0.1;                      //   GGX alpha of the sampling lobe (only the noise depends on it)
```

## 3. `src/Cli/SceneFileMaterials.cpp` — ParseMaterial

Add just before the `type == "emission"` branch:

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

Also add `measured` to the error message at the top (`"must be an object with a \"type\" (diffuse, conductor, ..."`)
and to the type list in the final `else Fail(...)`.

## 4. `src/Cli/Scenes/FileScene.cpp` — creating the material

Add the include:

```cpp
#include "Materials/MeasuredHpBRDF.hpp"
```

Add a case to `switch (m.type)` (next to `case MaterialSpec::Type::Conductor`):

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

`FileScene.cpp` is the only place that switches on `MaterialSpec::Type` (`grep -rn 'case MaterialSpec::Type'`).
If `/W4` warns about a missing enumerator, add `Measured` there.

## 5. One row in the material table of `docs/SCENE_FILE.md`

```
| `measured` | `file` (a PCA-compressed hpBRDF table `.hpbpca`, made with `python -m brdfkit export`), `alphaSample` (GGX alpha of the importance-sampling lobe; changes only the noise, not the image. Default 0.1) |
```

## 6. Checking it works

```bat
GoLD_rayt.exe scenes\hpbrdf_aluminum.json --spp 16 --width 400
GoLD_rayt.exe scenes\hpbrdf_aluminum.json --stokes --polarizer 0,90
```

- Loading prints `[SceneFile] material "aluminum": measured hpBRDF ...`. The RGB build takes a few seconds to bake the RGB table (OpenMP parallel).
- In the `--stokes` DoLP image, the degree of polarization should be high at the rim of the aluminium sphere (oblique incidence) and near 0 at the centre;
  then the Mueller matrices are passed correctly (same trend as the gold sphere).
- With polarizers at 0° and 90°, the top/bottom of the sphere (vertical plane of incidence) and the left/right swap brightness.

## 7. Spectral build

Rebuilding with `RAYT_ENABLE_SPECTRAL 1` in `Core/RenderMode.hpp` makes eval decode the table at each of the path's 4 wavelengths.
Of GoLD_rayt's 360–830 nm range, 360–414 nm lies outside the table and is extended with the 414 nm value (aluminium is nearly flat there, so the effect is small).
