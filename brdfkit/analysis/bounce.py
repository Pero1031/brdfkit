"""
Where three numbers per colour stop being enough: inter-reflection between the measured plastics, spectral against RGB.

A renderer that carries RGB multiplies a light's RGB by a surface's RGB at every bounce. The spectrum is multiplied per
wavelength, and the colour is taken at the end; the two agree only when the spectra are smooth. Dyed plastics have steep
absorption edges, so after a few bounces between coloured walls, or under a light with narrow bands, the RGB product drifts.
This module takes the reflectance spectra of the plastics from the reductions of `brdfkit quirks` and computes, per number of
bounces and light, the colour that the spectrum gives and the colour that the RGB product gives:

    reflectance   rho(lambda) of a sample: its intensity beside the lobe relative to a white sample's, times that white's albedo
    lights        spectra of a daylight (6500 K), an incandescent (2856 K) and a three-band screen
    compare       spectral and RGB colour after reflections off a sequence of surfaces: Delta E 2000 and the luminance ratio

RGB here is linear sRGB (negative components of a surface clamped at 0, as renderers do) or linear Rec. 2020. docs/BOUNCE.md.
"""
import numpy as np
from ..color import cmf, XYZ_TO_SRGB
from . import eyes as EY

XYZ_TO_REC2020 = np.array([[1.7166512, -0.3556708, -0.2533663],
                           [-0.6666844, 1.6164812, 0.0157685],
                           [0.0176399, -0.0427706, 0.9421031]])
SPACES = {'sRGB': XYZ_TO_SRGB, 'Rec. 2020': XYZ_TO_REC2020}


def reflectance(red, white, reference, white_albedo=0.65, theta_d_deg=(30.0, 45.0), lam_max=780.0):
    """(channels, rho) of a sample: floor-removed intensity beside the lobe relative to the white sample's, mean over theta_d,
    times `white_albedo`, clipped to [0, 1]; visible channels only"""
    s, w = EY.signals(red, reference, 'mid', theta_d_deg), EY.signals(white, reference, 'mid', theta_d_deg)
    lam = s['channels']; use = lam <= lam_max
    rho = white_albedo * (s['R'] / w['R']).mean(0)
    return lam[use], np.clip(rho[use], 0.0, 1.0)


def blackbody(lam, T):
    x = np.asarray(lam, float) * 1e-9
    return 1.0 / (x ** 5 * (np.exp(1.438777e-2 / (x * T)) - 1))


def lights(lam):
    """name -> spectrum (arbitrary scale): daylight as a 6500 K blackbody, incandescent 2856 K (CIE A), and a screen with three
    bands (Gaussian primaries at 455, 535, 625 nm, widths 10, 17, 12 nm, balanced to the daylight's colour)"""
    lam = np.asarray(lam, float)
    g = lambda mu, sg: np.exp(-0.5 * ((lam - mu) / sg) ** 2)                     # noqa: E731
    prim = np.stack([g(455, 10), g(535, 17), g(625, 12)], 1)                    # [C, 3]
    target = (blackbody(lam, 6500.0) @ cmf(lam))
    wts = np.linalg.solve((prim.T @ cmf(lam)).T, target)                         # primaries' weights giving the daylight's XYZ
    return {'daylight 6500 K': blackbody(lam, 6500.0), 'incandescent 2856 K': blackbody(lam, 2856.0), 'three-band screen': prim @ wts}


def _lab(xyz, white):
    t = np.asarray(xyz, float) / white
    f = np.where(t > (6 / 29) ** 3, np.cbrt(np.maximum(t, 0)), t / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def delta_e00(lab1, lab2):
    """CIEDE2000 colour difference"""
    L1, a1, b1 = lab1; L2, a2, b2 = lab2
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2); Cm = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cm ** 7 / (Cm ** 7 + 25.0 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1 = np.degrees(np.arctan2(b1, a1p)) % 360; h2 = np.degrees(np.arctan2(b2, a2p)) % 360
    dL, dC = L2 - L1, C2p - C1p
    dh = h2 - h1
    if C1p * C2p == 0:
        dh = 0.0
    elif dh > 180:
        dh -= 360
    elif dh < -180:
        dh += 360
    dH = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dh) / 2)
    Lm, Cpm = (L1 + L2) / 2, (C1p + C2p) / 2
    if C1p * C2p == 0:
        hm = h1 + h2
    elif abs(h1 - h2) <= 180:
        hm = (h1 + h2) / 2
    else:
        hm = (h1 + h2 + 360) / 2 if h1 + h2 < 360 else (h1 + h2 - 360) / 2
    T = 1 - 0.17 * np.cos(np.radians(hm - 30)) + 0.24 * np.cos(np.radians(2 * hm)) + 0.32 * np.cos(np.radians(3 * hm + 6)) \
        - 0.20 * np.cos(np.radians(4 * hm - 63))
    SL = 1 + 0.015 * (Lm - 50) ** 2 / np.sqrt(20 + (Lm - 50) ** 2); SC = 1 + 0.045 * Cpm; SH = 1 + 0.015 * Cpm * T
    RT = -2 * np.sqrt(Cpm ** 7 / (Cpm ** 7 + 25.0 ** 7)) * np.sin(np.radians(60 * np.exp(-((hm - 275) / 25) ** 2)))
    return float(np.sqrt((dL / SL) ** 2 + (dC / SC) ** 2 + (dH / SH) ** 2 + RT * (dC / SC) * (dH / SH)))


def compare(lam, light, rhos, space='sRGB'):
    """light reflected once off each surface of `rhos` in turn (a list of reflectance spectra): the spectral result against the
    RGB product. Both are scaled by the same factor, the one that brings the spectral result to the luminance of a mid grey
    (Y = 0.18 of the light), so that the colour difference is one of chromaticity and relative brightness, as after an exposure
    set for the patch. -> dict(dE00, Y_ratio (RGB / spectral luminance), xy_spectral, xy_rgb, clamped)"""
    A = SPACES[space]; Ai = np.linalg.inv(A); C = cmf(lam)
    white = light @ C
    spec = light.copy(); rgb = A @ white; clamped = False
    for rho in rhos:
        spec = spec * rho
        m = A @ (rho @ C) / (A @ C.sum(0))                    # the surface's RGB: its colour under an equal-energy light, white = 1
        clamped = clamped or bool((m < 0).any())
        rgb = rgb * np.maximum(m, 0.0)
    xyz_s, xyz_r = spec @ C, Ai @ rgb
    k = 0.18 * white[1] / xyz_s[1]
    xy = lambda v: tuple(v[:2] / v.sum())                     # noqa: E731
    return dict(dE00=delta_e00(_lab(k * xyz_s, white), _lab(k * xyz_r, white)), Y_ratio=float(xyz_r[1] / xyz_s[1]),
                xy_spectral=xy(xyz_s), xy_rgb=xy(xyz_r), clamped=clamped, level=float(xyz_s[1] / white[1]))


def report(names, lam, rhos, bounces=(1, 2, 3, 4, 6), say=print):
    L = lights(lam)
    say('Spectral against RGB after n reflections off the same plastic: Delta E 2000 (luminance of the RGB result / spectral)')
    say('each result is exposed to a mid grey before comparing; reflectance = intensity beside the lobe relative to the white sample x its albedo')
    for space in SPACES:
        for ln, light in L.items():
            say(f'\n{space}, {ln}' + ' ' * max(1, 34 - len(space) - len(ln)) + ''.join(f'{f"n = {n}":>16s}' for n in bounces))
            for name, rho in zip(names, rhos):
                cells = []
                for n in bounces:
                    r = compare(lam, light, [rho] * n, space)
                    cells.append(f'{r["dE00"]:6.1f} ({r["Y_ratio"]:4.2f})' + ('*' if r['clamped'] else ' '))
                say(f'  {name[:30]:32s}' + ''.join(f'{c:>16s}' for c in cells))
    say('\n* a component of the surface\'s RGB was negative and clamped at 0 (outside the gamut)')
    say('\none reflection off each of two different plastics (sRGB): Delta E 2000 (luminance ratio)')
    for ln, light in L.items():
        say(f'  {ln}')
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                r = compare(lam, light, [rhos[i], rhos[j]], 'sRGB')
                if r['level'] > 1e-4:
                    say(f'    {names[i][:22]:22s} then {names[j][:22]:22s} {r["dE00"]:6.1f} ({r["Y_ratio"]:4.2f})   light left {100 * r["level"]:.2f} %')
