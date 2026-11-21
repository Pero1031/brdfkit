# Corrected tables for Mitsuba 3

`brdfkit to-mitsuba` writes any table brdfkit can read (a measured hpBRDF table, or a factor model that was calibrated,
edited, swapped or roughened) as a `tensor_file` that the `measured_polarized` plugin of Mitsuba 3 loads, with the
corrections of [QUIRKS.md](QUIRKS.md) applied on the way. The 13 GB are streamed, one φ_d slice at a time.

```bash
# once: the floor and the factor g(theta_d), both from black glass
python -m brdfkit floor     data/kaist/hpbrdf/Black_glass.hpbrdf --out out/floor/Black_glass_floor.npz
python -m brdfkit factor    data/kaist/hpbrdf/Black_glass.hpbrdf --out out/clean/Black_glass_xc.npz \
                            --floor out/floor/Black_glass_floor.npz --exchange --cosine light
python -m brdfkit calibrate out/clean/Black_glass_xc.npz --out out/clean/black_glass_xc_cal.npz

# per material
python -m brdfkit to-mitsuba data/kaist/hpbrdf/White_smooth_plastic.hpbrdf out/clean/White_smooth_plastic.pbsdf \
       --floor out/floor/Black_glass_floor.npz --exchange --cosine light --calib out/clean/black_glass_xc_cal.npz

# a factor model (already cleaned when it was built): as it is
python -m brdfkit to-mitsuba out/gloss/Aluminum_a0.21.npz out/clean/Aluminum_a0.21.pbsdf --phi-step 4 --band-step 2
```

```xml
<bsdf type="measured_polarized">
    <string name="filename" value="out/clean/White_smooth_plastic.pbsdf"/>
</bsdf>
```

Implementation in `brdfkit/analysis/clean.py` (the corrections, as views of a table) and `brdfkit/io/hpbrdf.py`
(`write_table`).

## The corrections

What an hpBRDF table holds, as far as its own content shows:

    table(φ_d, θ_d, θ_h) = [ BRDF × cos(light) × g(θ_d) ] filed at φ_d + 180°  +  an unpolarized floor

| Option | What it does | Evidence | Where it holds |
|---|---|---|---|
| `--floor FLOOR.npz` | Subtracts the per-wavelength unpolarized offset measured on black glass from m00 | [FLOOR.md](FLOOR.md) | The floor differs between samples (0.6–0.9 of the glass's on the plastics, [UMOV.md](UMOV.md)), so on a dark sample the glass's floor can take too much; m00 is never taken below 2 % |
| `--exchange` | Shifts φ_d by 180°: the two directions trade places | [SLOTS.md](SLOTS.md) | Decided on the two samples with a body that were checked; no effect on a specular lobe |
| `--cosine light` | Divides by the cosine of the light's direction (clamped at `--min-cos`, 0.1) | Quirk 10 | Smooth samples' lobes, white plastics' bodies. Too much for the red rough plastic's body and the rough plastics' lobes |
| `--calib CAL.npz` | Divides by K(θ_d) of a calibration file | Quirks 5–6, [CALIBRATION.md](CALIBRATION.md) | Up to θ_d ≈ 50°–55° |

They are applied in that order, and each is optional. `--cosine` names the direction as labelled in the table that comes
out, so after `--exchange` it is `light`. A calibration file stores what the glass's lobe still lacks after the corrections
its model was built with, so it has to be used with the same `--exchange` and `--cosine`; another combination is refused.
With `--cosine light` the K in the file is g(θ_d) (0.50 at 20°, 0.36 at 45°, 0.29 at 60°–70°), which applies to lobes and
bodies; without it K also holds the cosine of the specular geometry (0.47 to 0.08) and is right for specular lobes only.

`--phi-step N` keeps every N-th φ_d bin and `--bands LO,HI` / `--band-step N` thin the wavelengths: the full table is
13 GB, `--phi-step 4 --band-step 2` gives 1.6 GB. θ_d and θ_h are kept, because the lobes are only a few bins wide.

## Checks

**The file.** Mitsuba 3.7.1 (`scalar_spectral_polarized`, the stock `measured_polarized`) loads the written files. At bin
centres its `eval` agrees with the table, rotated from the table's frames into Mitsuba's Stokes bases with
`brdfkit.geometry` and multiplied by the light's cosine as the plugin does, to 5 × 10⁻⁶ of m00 in all 16 elements on a
synthetic table that is not symmetric in light and camera, as stored and exchanged, and to 10⁻³ on the corrected white
plastic. Bins below the horizon are NaN, as in the original files; no evaluation above the horizon returned a NaN
(8000 random directions, a quarter of them within 4° of grazing), and a sphere under constant light renders without one.

**The corrections on real tables.** The white smooth plastic, corrected:

| | As stored | Corrected |
|---|---|---|
| `brdfkit slots` verdict | exchanged | as labelled |
| Orientation of the body's polarization against the prediction, as labelled | −0.13 | +0.98 |
| m00 of the body ∝ cos^a("light") cos^b("camera"), as labelled | a = −0.04, b = +0.78 | a = −0.16, b = −0.05 |
| Directional albedo at 15°–55° (luminance) | 0.22–0.25 | 0.68–0.75 |

Aluminium, whose lobe was not used for the calibration: directional albedo / Fresnel (Rakić, 20 nm film) at incidence
25°, 35°, 45°, 55°, 65°, 75° is 0.43, 0.37, 0.27, 0.19, 0.13, 0.10 with only the floor removed and
**1.08, 1.09, 1.01, 0.91, 0.84, 0.87** corrected (`calibrate --check`). The route without the cosine
([CALIBRATION.md](CALIBRATION.md)) gives 1.01, 1.05, 1.00, 0.92, 0.85, 0.86: for a specular lobe the two are the same
correction within a few percent.

## Limits

- The scale is restored to about 50°–55° of θ_d. Beyond that the corrected tables are too bright by a factor that differs
  between materials (the white plastics' albedo passes 1 at an incidence of 75°); for stimuli or comparisons keep light and
  camera less than about 110° apart.
- The cosine of the light is clamped at 0.1 (84°), so the last degrees before grazing incidence are under-corrected rather
  than amplified.
- Nothing is done about the lobe's width (the instrument's angular response, 1.1° to 3°), the filled-in ends of θ_d, the
  irregular specular spectrum, or the metals' first row and column (quirks 1–3, 7, 9). For a metal's polarization use a
  factor model with N from its n, k (`brdfkit swap --nk`) and export that.
- A factor model keeps the polarization as a function of θ_d and wavelength only. That is the specular lobe; the
  polarization of a body, which depends on the two directions separately, needs the table itself.
- hpBRDF has no stated licence. The command writes a file on your disk; do not redistribute corrected tables.
