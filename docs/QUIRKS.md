# Quirks of the hpBRDF tables

The hpBRDF release (Moon et al., SIGGRAPH Asia 2025, arXiv 2509.13779) gives 14 materials as tables of Mueller matrices:
361 φ_d × 91 θ_d × 91 θ_h bins, 68 wavelengths from 414 to 950 nm. Anyone who renders them, fits a model to them or compares
them with other data should know where they depart from the reflectance of a physical surface. This page collects those
departures for eight materials: aluminium, SUJ2 steel, fake gold, black glass, and four plastics (white smooth, white rough,
red rough, black rough). All of them were measured the same way, by one command. Most of the departures appear on every material and come from the
measurement or from how the tables were built; two belong to the samples.

```bash
python -m brdfkit quirks data/kaist/hpbrdf/Black_glass.hpbrdf --out out/quirks/Black_glass.npz --dielectric
python -m brdfkit quirks data/kaist/hpbrdf/Aluminum.hpbrdf --out out/quirks/Aluminum.npz --nk data/kaist/refidx/Al_Rakic.yml --film-nm 20
python -m brdfkit quirks out/quirks/Aluminum.npz --out out/quirks/Aluminum.npz --nk ...    # report again from the reduction, in a second
```

One pass reads the 13 GB file (2 to 6 minutes) and keeps a 12 MB reduction: the luminance of every bin, and Mueller
matrices summed over φ_d at the specular core (θ_h ≤ 0.55°), beside the lobe (θ_h 15°–25°, the body of a diffuse material)
and far from it (θ_h ≥ 40°). `--dielectric` fits the index from the specular polarization; `--nk` takes a metal's n, k.
Implementation in `brdfkit/analysis/quirks.py`.

The instrument lights the sample at 40° to 160° from the camera in 5° steps, so the measured θ_d (half that angle) runs from
about 20° to 80° in 2.5° steps. The checks below use θ_d 20°–79°.

## Summary

| # | Quirk | Size on the eight materials | Origin | What to do |
|---|---|---|---|---|
| 1 | Both ends of θ_d are filled in | The specular value levels off below 7°–14° and above 79°–89° (less than 0.1 % change per bin) | How the tables were built (Gaussian inpainting, stated in the paper) | Do not fit or compare there; `brdfkit factor` copies the edge values |
| 2 | The measured range is smoothed too | Bin-to-bin roughness of the specular value along θ_d 0.003–0.015 (MERL: about 0.02) | Writing each pixel into neighbouring bins, and the same Gaussian fill | Nothing; keep it in mind for fine angular structure |
| 3 | The peak is flat and the lobe widens with θ_d | Value at θ_h 0.28° / at 0: 0.97–0.99 on the smooth samples. These all have the same half width: 1.1° at θ_d 30°, 3.0° at 75° | The instrument's angular response | A lobe narrower than this is not resolved; `brdfkit gloss` builds on it |
| 4 | An unpolarized floor far from the lobe | 0.007–0.010 at 550 nm on dark samples (it differs by about ±15 % between samples); 4–5 × that at 414 nm, 2 × at 950 nm; rises toward grazing | The instrument (an additive offset) | `brdfkit floor` ([FLOOR.md](FLOOR.md)) |
| 5 | The specular lobe carries far too little energy | Lobe energy / Fresnel 0.22–0.36 at θ_d 30°, 0.03–0.10 at 75°, the same on metals, glass and plastic; a rough black plastic holds about a tenth | Unknown; specific to the surface reflection | `brdfkit calibrate` for specular materials, up to θ_d ≈ 50° ([CALIBRATION.md](CALIBRATION.md)) |
| 6 | The diffuse body is dark, by another factor | White plastics, smooth and rough: directional albedo 0.21–0.26 from the table, at every incidence | Unknown | No correction in the data set; use the body in relative terms |
| 7 | The specular spectrum varies irregularly along θ_d | 0.7 %–3.7 % per bin; far from the lobe only 0.1 %–0.3 % | Specific to the specular peak | Smooth S along θ_d or use a physical model ([FACTOR.md](FACTOR.md)) |
| 8 | One-sided φ_d arms | Arm / its mirror image at θ_d 25°: SUJ2 7.0, aluminium 3.4, fake gold 3.3, black glass 2.3, red plastic 1.4; white and black plastics 1.1 (none) | The samples (polishing marks) | `--sym` in `factor` and `edit` ([LATENT_EDIT.md](LATENT_EDIT.md)) |
| 9 | Metals' polarization disagrees with Fresnel; dielectrics' agrees | Diattenuation at θ_d 60°: aluminium 0.28 (Fresnel 0.06), SUJ2 0.02 (iron 0.44); five dielectrics within 0.07 (0.1 at 45°) | Metals only | Take N of a metal from its n, k (`brdfkit swap --nk`) |

Smaller items: reciprocity (φ_d 0° against 180°) holds within 20 % at θ_d 25°–45° and is off by up to 57 % (either way) at
θ_d ≥ 65°; the lobe energy fluctuates by 3 %–11 % from one θ_d bin to the next around a smooth trend.

## Numbers per material

At four values of θ_d (550 nm for the polarization):

| | θ_d | Black glass | Aluminium | SUJ2 | Fake gold | White smooth plastic | White rough plastic | Red rough plastic | Black rough plastic |
|---|---|---|---|---|---|---|---|---|---|
| Lobe energy / Fresnel | 30° | 0.356 | 0.354 | 0.276 | – | 0.333 | 0.218 | 0.227 | (0.085) |
| | 45° | 0.227 | 0.200 | 0.182 | – | 0.193 | 0.139 | 0.139 | (0.057) |
| | 60° | 0.132 | 0.088 | 0.124 | – | 0.109 | 0.090 | 0.073 | (0.039) |
| | 75° | 0.075 | 0.062 | 0.103 | – | 0.062 | 0.049 | 0.033 | (0.016) |
| Half width of the lobe [°] | 30° | 1.15 | 1.76 | 1.53 | 1.11 | 1.12 | 4.02 | 4.60 | 8.32 |
| | 75° | 3.06 | 4.26 | 3.61 | 2.81 | 3.00 | 4.88 | 6.43 | 8.82 |
| m01 (row) / m10 (column) / Fresnel | 45° | 0.73 / 0.73 / 0.82 | 0.16 / 0.07 / 0.03 | −0.05 / −0.15 / 0.23 | 0.18 / 0.09 / – | 0.73 / 0.70 / 0.79 | 0.77 / 0.75 / 0.82 | 0.82 / 0.80 / 0.83 | 0.79 / 0.80 / 0.85 |
| | 60° | 0.95 / 0.96 / 0.98 | 0.28 / 0.23 / 0.06 | 0.02 / −0.13 / 0.44 | 0.19 / 0.14 / – | 0.96 / 0.96 / 1.00 | 0.94 / 0.88 / 0.99 | 0.96 / 0.96 / 0.98 | 0.90 / 0.96 / 0.97 |
| Index from the polarization | | 1.533 | – | – | – | 1.624 | 1.555 | 1.505 | 1.457 |

Fresnel is the fitted index for the dielectrics, Rakić aluminium under a 20 nm oxide film, and iron (Johnson & Christy) for
SUJ2; fake gold has no known index. The lobe energy is integrated over θ_h < 10° above the median beside the lobe
(θ_h 20°–30°), which takes out the body of the plastics. The rough plastics' wider lobes reach past that baseline, so their
ratios read lower; for the black rough plastic (half width 8°–9°) the measure does not hold, hence the brackets (see 5–6).

## 1–2. Filled-in and smoothed angles

The paper fills the bins the measurement did not reach (mostly low θ_d) by a Gaussian convolution over (φ_d, θ_d, θ_h), per
Mueller element and wavelength. In the tables the specular value stops changing toward both ends of θ_d. The edge is soft,
because the fill blends into the measured range. On black glass the value changes by less than 10⁻⁴ per bin below 8°, 10⁻³
below 10° and 10⁻² below 14°, and likewise above 89°, 86° and 82°; the other materials are within a few degrees of that.
Anything that compares hpBRDF with another data set at normal incidence or at grazing compares it with the fill. An earlier
reading of the aluminium's polarization at normal incidence was withdrawn for this reason ([LU_CHIPMAN.md](LU_CHIPMAN.md)).

Inside the measured range the tables are also smoother than a direct measurement: the specular value's bin-to-bin roughness
along θ_d (median |second difference| of its log) is 0.003–0.015 against about 0.02 for the MERL tables, which are averages
of raw samples. In the filled-in range it falls to 10⁻⁵–10⁻⁴.

## 3. Flat peak, and the instrument's angular response

On the smooth samples the value at θ_h 0.28° is 97 %–99 % of the value at θ_h 0, even on polished black glass, where the
reflection is a mirror's (the rough plastics' lobes are flat there anyway). One camera pixel on the 410-pixel image of the sample sphere spans at least 0.28° of surface
normal, and each pixel is also written into the neighbouring bins, so the first eight θ_h bins (0°–0.54°) come from the same
pixels.

The half width of the lobe is the same on the three smooth samples, black glass, white smooth plastic and fake gold: 1.1° at
θ_d 30° (within 0.05°) and 2.8°–3.1° at 75°. A surface's lobe does not widen like that; the instrument's angular response does.
Measured lobes are that response convolved with the surface's; [GLOSS.md](GLOSS.md) uses this.

## 4. The floor

Far from the lobe (θ_h ≥ 40°), the dark samples read 0.007–0.010 at 550 nm: black rough plastic 0.0068, black glass 0.0084,
SUJ2 0.0093, fake gold 0.0095 (the per-bin median used in [FLOOR.md](FLOOR.md) gives 0.006 on black glass). So the floor is
nearly, not exactly, the same from sample to sample (about ±15 %), and subtracting the glass's floor from another dark sample
can take too much. It carries no polarization (every |m_ij / m00| below 0.03 on
black glass), and has a U-shaped spectrum: 4.2–4.9 times the 550 nm value at 414 nm and about 2 times at 950 nm on all four dark samples. It rises
toward grazing (1.5 times at θ_d 75° on black glass). On the red rough plastic, which absorbs blue, the floor is what the
table holds at 414 nm (414 / 550 nm = 4.3 far from the lobe). In a render the floor turns black glass navy blue
([CALIBRATION.md](CALIBRATION.md)).

## 5–6. The radiometric scale: the lobe, and the body

**The lobe.** Divided by the Fresnel reflectance of the sample's own index, the energy of the specular lobe is about a third
at θ_d 30° and a tenth or less at 75°. The ratio is the same on black glass, aluminium, white smooth plastic and (within
about 20 %) SUJ2, materials whose reflectance differs by a factor of 20 and whose lobes differ in width; the rough white and
red plastics read a little lower because part of their wider lobe falls into the baseline. The black rough plastic, whose lobe
is 8°–9° wide, is beyond this measure, but its whole table holds a directional albedo of only 0.016–0.023 at 25°–75°,
about what the floor alone gives, where Fresnel of its index (1.46) gives 0.035–0.24. Its lobe above the minimum beside it
carries roughly a tenth of Fresnel at 45°. So a rough surface reflection is short of energy as well. So it belongs to the
tables, not to the samples. Up to θ_d ≈ 50° one correction K(θ_d), measured on black glass, serves all of them within
10 %–15 %; beyond 55° each material keeps its own error of 20 %–35 % ([CALIBRATION.md](CALIBRATION.md)).

**The body.** The white plastics have a diffuse body under their specular lobe. From the table, by Monte Carlo over the
hemisphere, the directional albedo (luminance) of the smooth one is (the rough one is within 0.03 of it, 0.21–0.25):

| Incidence | 15° | 25° | 35° | 45° | 55° | 65° | 75° |
|---|---|---|---|---|---|---|---|
| As delivered | 0.22 | 0.23 | 0.23 | 0.24 | 0.25 | 0.26 | 0.26 |
| Body only (θ_h ≥ 10°) | 0.17 | 0.18 | 0.20 | 0.22 | 0.23 | 0.25 | 0.24 |
| Divided by the glass's K(θ_d) | 0.54 | 0.57 | 0.62 | 0.70 | 0.81 | 0.96 | 1.10 |

An ordinary white plastic reflects 0.7–0.9, so the body is probably 3 to 4 times too dark, the same on both white plastics. But it is dark by a nearly
constant factor: dividing the whole table by the lobe's K(θ_d) over-corrects it, up to an albedo above 1 at 75°. At a fixed
geometry (light and camera both 45° from the normal, θ_d set by the azimuth between them) the body changes by −17 % from
θ_d 21° to 43°, where K changes by −43 %. So the strong θ_d dependence of K is specific to the specular peak, and
`brdfkit calibrate` must not be applied to a material with a diffuse body. The data set has no diffuse reference (a
Spectralon-like sample) that would calibrate the body; the paper's supplementary document may describe the radiometric
calibration.

## 7. Spectrum along θ_d

The specular spectrum (normalized per θ_d) changes irregularly from one θ_d bin to the next: by 0.7 % (red plastic) to 3.7 %
(fake gold), as the RMS second difference of its log, and on black glass, which reflects almost the same at every wavelength,
the specular spectrum is tilted (414 nm 1.27, 950 nm 0.68 relative to 550 nm at θ_d 30°). The spectrum far from the lobe
(the body, or the floor) is smooth along θ_d on every material: 0.1 %–0.3 %. The two do not vary together (correlation
−0.2 to +0.2). So the spectral calibration per illumination angle is consistent, and the irregularity comes with the bright
specular peak itself. The factor model keeps S(θ_d, λ) as measured; smoothing it or replacing it by a physical model removes
the irregularity ([FACTOR.md](FACTOR.md), [CALIBRATION.md](CALIBRATION.md)).

## 8. φ_d arms and reciprocity

An isotropic BRDF is the same at φ_d and −φ_d. Most samples have a one-sided "arm" on the shoulder of the lobe (θ_h 6°–15°):
SUJ2 at φ_d −52° with 7.0 times the value at +52° (and 7.6 times at θ_d 45°), aluminium at −24° (3.4 times), fake gold at
−147° (3.3), black glass at −121° (2.3), red rough plastic at −117° (1.4). The white smooth, white rough and black rough plastics have
none (1.1). The direction differs from sample to sample, and moulded plastics have none, so the arms are the samples' own anisotropy
(polishing marks on the metals, scratches on the glass), not the instrument's. They show as streaks from the highlight once a
lobe is widened ([LATENT_EDIT.md](LATENT_EDIT.md)). Averaging the four equivalent angles φ, −φ, φ+180°, 180°−φ removes them
(`--sym median`).

Reciprocity asks the same of φ_d and φ_d + 180°. At θ_d 25°–45° the two agree within 20 % (θ_h 1°–10°). At θ_d ≥ 65°
φ_d 0° reads 1.03–1.57 times φ_d 180° on seven materials and 0.72 times on the black rough plastic.

## 9. Polarization at the specular

For an ideal reflection, the diattenuation (m01, a row of the Mueller matrix in the (s, p) bases) equals the polarizance
(m10, a column), and both follow Fresnel. The five dielectrics do: m01 and m10 agree within 0.01–0.07 at θ_d 45°–60°,
reach 0.90–0.96 at Brewster's angle (Fresnel 0.97–1.00), and give plausible indices (black glass 1.53, plastics 1.46–1.62). They read up to 0.1 below
Fresnel at θ_d 45°, and the black glass reads far lower below 35° (0.19 against 0.39 at 30°): the instrument error that
[CALIBRATION.md](CALIBRATION.md) measured on the glass.

The metals do not: aluminium reads 0.28 at θ_d 60°, five times its Fresnel value under a film (0.06), and m10 differs from
m01; SUJ2 reads 0.02 where iron's Fresnel gives 0.44, and its m10 is negative over θ_d 30°–65°. Fake gold has no index to
compare with, but its m01 and m10 differ in the same way. So the first row and column of the metals' Mueller matrices are not
usable as measured. The polarization a metal gives the camera under unpolarized light, and what it does to polarized light
(a polarized sky, a second reflection), should come from its n, k. The conclusion of [LU_CHIPMAN.md](LU_CHIPMAN.md) that the
aluminium's diattenuation exceeds theory at θ_d 60°–79° is part of this problem, not a property of the sample.

## What to do with the tables

- Use θ_d 20°–80° only; outside it the tables are a fill.
- Remove the floor before looking at anything dark: dark materials, the far field, the blue and infrared ends.
- Do not read absolute reflectance from the tables. For specular materials, `calibrate` restores it up to θ_d ≈ 50°; for
  materials with a diffuse body there is no correction in the data set.
- Expect lobes no narrower than about 1° (θ_d 30°) to 3° (75°), whatever the surface.
- Symmetrize over φ_d before widening or editing a lobe.
- Smooth the specular spectrum along θ_d, or model it.
- Trust the polarization of dielectrics (to about 0.1, better near Brewster's angle), not that of metals.

## Limits

Eight of the 14 materials were checked: three metals, black glass and four plastics. The two silicones, the billiard ball and
the plum, silver and yellow rough plastics are not; the four plastics checked cover smooth against rough and a bright against
a dark body. The causes of quirks 5, 6, 7 and 9 are not known; the paper's supplementary document, not
available from the authors' page, may describe the calibration. The white plastic's true albedo is not known either: "3 to 4
times too dark" assumes an ordinary white plastic.
