# The table as an ellipsometer: what is the fake gold made of?

One hpBRDF sample, "Fake gold", has no known optical constants. A table of Mueller matrices holds, at its specular geometry,
what an ellipsometer measures, at every angle of incidence and 68 wavelengths. This page asks how far that goes: it checks
the method on two metals whose constants are known, then applies it to the fake gold. The answer is that the known metals
are recognized and the fake gold is not identified: it is a gold-coloured conductor that fits brass, titanium nitride, gold
and copper about equally, and none of them well.

```bash
python -m brdfkit ellipso out/quirks/Fake_gold.npz --candidates data/kaist/refidx/candidates/*.yml \
       --factor out/factor/Fake_gold_s0.npz --calib out/calib/black_glass_cal.npz
```

It works on the reduction of `brdfkit quirks` (the specular core) and, for the colour and the level, on the factor model and
the black-glass calibration. Candidates are n, k tables from refractiveindex.info (CC0). Implementation in
`brdfkit/analysis/ellipso.py`.

## What can be read

A reflection off a smooth surface is ρ = r_p / r_s = tan Ψ · e^{iΔ}. In the normalized Mueller matrix, in (s, p) bases:

    m01 = m10 = −cos 2Ψ,    (m22 + m33) / 2 = sin 2Ψ cos Δ,    (m23 − m32) / 2 = sin 2Ψ sin Δ

Ψ sits in the first row and column, which the hpBRDF tables get wrong on metals ([QUIRKS.md](QUIRKS.md), 9). **Δ does not
need them**: it is the angle of the last two terms, and depolarization, which shrinks both, leaves it alone. So the phase is
the part of a metal's polarization that can be used. It falls from 180° at normal incidence through 90° at the principal
angle, and the curve depends on n and k.

Two more things are available at the specular, both less exact: the spectrum of the reflection (S of the factor model,
good to 10 %–25 %), and its level, the lobe's energy divided by the black glass's K (good to 10 %–15 % up to θ_d 50°).

## Check 1: is the phase right?

Aluminium against Rakić's constants, bare and under a 20 nm film of index 1.65 (the film found in
[LU_CHIPMAN.md](LU_CHIPMAN.md)); Δ in degrees:

| θ_d | 25° | 35° | 45° | 55° | 65° | 75° |
|---|---|---|---|---|---|---|
| Measured, 462 nm | 173.1 | 164.8 | 153.9 | 137.2 | 115.6 | 74.9 |
| Rakić + 20 nm | 172.0 | 163.9 | 152.2 | 135.7 | 112.0 | 77.4 |
| Rakić, bare | 176.0 | 171.9 | 165.9 | 156.9 | 142.6 | 115.7 |
| Measured, 550 nm | 174.0 | 166.3 | 156.6 | 137.5 | 121.8 | 85.6 |
| Rakić + 20 nm | 173.3 | 166.5 | 156.5 | 142.1 | 120.8 | 87.0 |

The phase follows the filmed aluminium within 4° at 462 nm; at longer wavelengths within 5° up to θ_d 50° and within 10°
beyond. The same matrices' m01 reads 0.2–0.3 against a Fresnel value of 0.05. Black glass switches from 180° to 0° between 55° and 60° (Brewster's angle
for 1.53: 56.8°). So Δ is usable; its scatter along θ_d is 2°–4°.

## Check 2: does the phase identify a metal? No

Ranking ten candidates by the RMS difference of Δ alone, each under a transparent film (index 1.5) of its best thickness,
aluminium comes out fifth:

| Candidate for the aluminium | Fe | Cr | Ti | Ni | Al | Ag | Cu |
|---|---|---|---|---|---|---|---|
| Best film [nm] | 17 | 20 | 19 | 16 | 31 | 11 | 5 |
| RMS of Δ [°] | 2.8 | 2.9 | 3.5 | 4.1 | 4.3 | 6.4 | 6.5 |

A thin film and the metal under it trade against each other: the phase gives roughly one number per wavelength (the
principal angle), and a film of unknown thickness uses it up. Fitting n and k per wavelength to Δ(θ_d) has the same
problem: k comes out plausible and n is not determined (0.7–1.7 on aluminium, pinned at 0 on the fake gold).

## Check 3: phase, colour and level together

The film changes the phase a lot and the reflectance little, so the colour and the level break the tie. "Combined" is the
RMS of the three misfits over their tolerances (3°, 0.15 and 0.15 in log units); about 1 is as good as the data allow.

| Aluminium (reflectance 0.78) | film [nm] | combined | Δ [°] | shape | level |
|---|---|---|---|---|---|
| **Al (Rakić)** | 29 | **1.11** | 4.0 | 0.14 | 0.15 |
| Ni | 13 | 1.42 | 3.7 | 0.21 | 0.24 |
| Ag | 8 | 1.53 | 5.6 | 0.15 | 0.24 |
| Ti, Cu, brass, Cr, Fe, Au, TiN | | 1.67–2.85 | | | |

| SUJ2 steel (reflectance 0.39) | film [nm] | combined | Δ [°] | shape | level |
|---|---|---|---|---|---|
| **Cr** | 58 | **1.77** | 7.9 | 0.20 | 0.12 |
| **Ti** | 61 | 1.84 | 7.7 | 0.27 | 0.08 |
| **Fe** | 40 | 1.92 | 9.1 | 0.20 | 0.04 |
| **Ni** | 61 | 2.00 | 8.0 | 0.31 | 0.13 |
| Cu, TiN, brass, Al, Au, Ag | | 2.90–3.92 | | | |

Aluminium is recognized as aluminium, by a modest margin (0.3). The bearing steel is recognized as one of the grey
transition metals, which cannot be told apart here, and clearly not as a noble metal or aluminium. Its phase fits poorly
(8°) and wants a film of 40–60 nm: its principal angle is 65°–68° (438–702 nm) where iron's is 76°–77°, so its surface is not bare steel
(an oxide or a preservative), which also explains why its measured polarization disagreed with iron in
[QUIRKS.md](QUIRKS.md).

## The fake gold

| Fake gold (reflectance 0.61) | film [nm] | combined | Δ [°] | shape | level | model reflectance |
|---|---|---|---|---|---|---|
| Brass (Cu70Zn30) | 12 | 1.59 | 3.5 | 0.28 | 0.25 | 0.77 |
| TiN | 4 | 1.71 | 4.3 | 0.29 | 0.26 | 0.48 |
| Ni | 27 | 1.72 | 4.6 | 0.37 | 0.10 | 0.58 |
| Au | 10 | 1.73 | 5.0 | 0.28 | 0.24 | 0.76 |
| Cu | 14 | 1.80 | 5.0 | 0.38 | 0.12 | 0.66 |
| Ti, Cr, Fe, Al, Ag | | 1.95–2.63 | | | | |

Four gold- or copper-coloured conductors, and nickel on the strength of its level alone, lie within 0.2 of each other,
less than the margin that separated aluminium from nickel. None reaches the 1.1 of a correct identification. What the data
do say:

- It is a conductor with a thin transparent film or none (4–14 nm for the gold-coloured candidates), not a dielectric
  coating that makes the colour by interference: the phase passes 90° like a metal's at every wavelength up to 650 nm.
- It is darker than gold or brass (0.61 against 0.77) and brighter than titanium nitride (0.48).
- Its principal angle rises more steeply with wavelength than any candidate's:

| Principal angle [°] | 438 nm | 462 | 502 | 550 | 598 | 646 |
|---|---|---|---|---|---|---|
| Fake gold | 61.7 | 64.1 | 69.6 | 71.1 | 73.5 | 73.9 |
| Brass + 12 nm | 64.4 | 64.5 | 66.4 | 68.9 | 71.0 | 72.7 |
| TiN + 4 nm | 64.2 | 64.7 | 66.0 | 68.7 | 71.0 | 72.6 |
| Au + 10 nm | 66.4 | 65.9 | 65.1 | 67.9 | 70.7 | 72.8 |
| Cu + 14 nm | 66.8 | 67.5 | 68.2 | 68.2 | 69.4 | 71.9 |

  12° from 438 to 646 nm, against 8° for brass and titanium nitride and 6° for gold. Its spectrum also keeps rising into
  the near infrared (1.30 at 646 nm, 1.64 at 798 nm relative to luminance 1) where gold and brass level off.

Both point to something that absorbs in the blue more than these metals do: another alloy (a brass with less zinc, a
bronze, a nitride of other composition), or a metal under a film that is not transparent, such as a tinted lacquer. Neither
is in the candidate set, and a tinted film is outside the model.

## Limits

- Ten candidates, one table of constants each. Alloys and nitrides vary with composition and deposition; the brass is one
  70/30 sample and the titanium nitride one film.
- The film is one transparent layer of fixed index 1.5. An absorbing or thick layer is not modelled.
- The tolerances are set from what the known samples show (the aluminium's spectrum is off by up to 27 % from a flat one),
  not from an error model, and the margins (0.3 for aluminium) are small. Two known metals are a thin validation.
- The level rests on the black-glass calibration, which holds to 10 %–15 % up to θ_d 50° ([CALIBRATION.md](CALIBRATION.md)).
- Ψ is not used at all. On an instrument whose first row and column were right, Ψ and Δ together would give n and k
  directly; here half of the ellipsometric information is missing.
