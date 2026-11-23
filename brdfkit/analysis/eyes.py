"""
Other eyes: what a near-infrared camera and a polarization camera see of the measured plastics, off the specular.

A camera that looks at a diffuse object under unpolarized light, away from the highlight, records per band the intensity R and
the linearly polarized part (q, u) = (m10, m20) of the first column of the Mueller matrix. Two things reach it: light reflected
by the surface (polarized along s, the pigment does not colour it) and light that went through the body (the pigment's colour,
weakly polarized along p by the refraction on the way out). So the polarization a camera sees is

    q(lambda) = q_surface - |q_body| B(lambda)

and it changes sign where the body gets bright enough: the angle of polarization turns by 90 degrees between a band where the
pigment absorbs and one where it does not. This module reads those quantities from the reductions of `brdfkit quirks`
(theta_h 15-25 deg, summed over phi_d, table frames: q > 0 is the surface's sign), with the instrument's unpolarized floor
taken from a dark reference (black glass: m00 - m11 in the same region), and answers for a set of materials:

    bands        R, q / R per camera band (blue, green, red, near infrared)
    crossing     the wavelength where q changes sign (NaN: it does not)
    contrast     R relative to a white reference, per band: what a near-infrared camera makes of the colours
    best_band    the band in which the weakest polarization among the materials is strongest

docs/EYES.md.
"""
import numpy as np

BANDS = {'blue': (430.0, 490.0), 'green': (510.0, 570.0), 'red': (600.0, 660.0), 'NIR': (800.0, 900.0)}
THETA_D = (30.0, 45.0, 60.0, 70.0)


def signals(red, reference, region='mid', theta_d_deg=THETA_D):
    """quirks reduction of a material and of a dark reference -> dict(channels, theta_d [deg], R [T, C] the intensity with the
    floor removed, q, u [T, C] the polarized part (m10, m20), floor [T, C])"""
    td = np.degrees(np.asarray(red['theta_d'], float))
    rows = [int(np.argmin(abs(td - t))) for t in theta_d_deg]
    with np.errstate(invalid='ignore', divide='ignore'):
        M = red[f'{region}_M'][rows] / red[f'{region}_n'][rows][:, None, None, None]
        G = reference[f'{region}_M'][rows] / reference[f'{region}_n'][rows][:, None, None, None]
    floor = G[..., 0, 0] - np.abs(G[..., 1, 1])
    return dict(channels=np.asarray(red['channels'], float), theta_d=td[rows], R=M[..., 0, 0] - floor, q=M[..., 1, 0], u=M[..., 2, 0],
                floor=floor)


def bands(sig, band_defs=BANDS):
    """per band: R [T, B], q [T, B], u [T, B] (means over the band's channels) and names"""
    lam = sig['channels']; out = dict(names=list(band_defs), theta_d=sig['theta_d'])
    for key in ('R', 'q', 'u', 'floor'):
        out[key] = np.stack([sig[key][:, (lam >= a) & (lam <= b)].mean(1) for a, b in band_defs.values()], 1)
    return out


def crossing(sig, lam_range=(430.0, 900.0), smooth=5, min_share=0.2):
    """per theta_d: wavelength [nm] of the first sign change of q (smoothed over `smooth` channels) inside lam_range, NaN if none;
    and the sign of q at the short end (+1: surface, -1: body). A change counts only when q reaches `min_share` of its largest
    magnitude on both sides of it: q that merely falls to zero (a bright body under a weak surface) is not a change of sign."""
    lam = sig['channels']; use = (lam >= lam_range[0]) & (lam <= lam_range[1])
    k = np.ones(smooth) / smooth
    out = np.full(len(sig['theta_d']), np.nan); first = np.zeros(len(sig['theta_d']))
    for t, q in enumerate(sig['q']):
        qs = np.convolve(q, k, mode='same')[use]; l = lam[use]
        first[t] = np.sign(qs[0]); lim = min_share * np.abs(qs).max()
        for i in np.nonzero(np.sign(qs[1:]) != np.sign(qs[:-1]))[0]:
            if np.abs(qs[:i + 1]).max() >= lim and (np.sign(qs[i + 1]) * qs[i + 1:]).max() >= lim:
                out[t] = l[i] + (l[i + 1] - l[i]) * qs[i] / (qs[i] - qs[i + 1]); break
    return out, first


def contrast(b, b_white):
    """R of a material relative to a white reference, [T, B]"""
    return b['R'] / b_white['R']


def best_band(bs, min_R=None):
    """bs: list of `bands` results. Per theta_d: the weakest |q| / R among the materials for each band [T, B], and the index of the
    band where that is largest. A material whose R is below `min_R` [B] in a band (at the floor's uncertainty) counts as fully
    polarized there: it is dark, so its surface is all there is to see."""
    p = np.stack([np.abs(b['q']) / np.maximum(b['R'], 1e-12) for b in bs])           # [M, T, B]
    if min_R is not None:
        dark = np.stack([b['R'] < min_R for b in bs])
        p = np.where(dark, 1.0, p)
    worst = np.minimum(p, 1.0).min(0)
    return worst, worst.argmax(1)


def report(names, sigs, white=None, say=print):
    bs = [bands(s) for s in sigs]; bn = bs[0]['names']; td = bs[0]['theta_d']
    say('What a camera sees off the specular (theta_h 15-25 deg, unpolarized light), floor removed')
    say('bands: ' + ', '.join(f'{k} {a:.0f}-{b:.0f} nm' for k, (a, b) in BANDS.items()))
    say('\nintensity R x 1000 per band' + (f', and relative to {names[white]}' if white is not None else ''))
    for t, th in enumerate(td):
        say(f'  theta_d {th:.0f}' + ' ' * 14 + ''.join(f'{k:>8s}' for k in bn) + ('      ' + ''.join(f'{k:>7s}' for k in bn) if white is not None else ''))
        for n, b in zip(names, bs):
            line = f'    {n[:22]:22s}' + ''.join(f'{1000 * v:8.1f}' for v in b['R'][t])
            if white is not None:
                line += '      ' + ''.join(f'{v:7.2f}' for v in contrast(b, bs[white])[t])
            say(line)
    say('\npolarization q / R [%] per band (+: the surface\'s direction, -: the body\'s, 90 deg from it); "dark": R below a third of the floor,')
    say('where the ratio is not reliable (the polarized part q x 1000 is given in brackets)')
    for t, th in enumerate(td):
        say(f'  theta_d {th:.0f}' + ' ' * 14 + ''.join(f'{k:>14s}' for k in bn) + '   sign change at')
        for n, b, s in zip(names, bs, sigs):
            cells = []
            for j in range(len(bn)):
                R, q = b['R'][t, j], b['q'][t, j]
                cells.append(f'dark ({1000 * q:5.2f})' if R < b['floor'][t, j] / 3 else f'{100 * q / R:6.1f} ({1000 * q:5.2f})')
            x, first = crossing(s)
            say(f'    {n[:22]:22s}' + ''.join(f'{c:>14s}' for c in cells)
                + (f'   {x[t]:.0f} nm ({"surface" if first[t] > 0 else "body"} below it)' if np.isfinite(x[t]) else '   -'))
    worst, pick = best_band(bs, min_R=None)
    say('\nweakest |q| / R among the materials [%] per band: the band a polarization camera should use is the one where this is largest')
    for t, th in enumerate(td):
        say(f'  theta_d {th:.0f}  ' + ''.join(f'{k} {100 * worst[t, j]:.1f}   ' for j, k in enumerate(bn)) + f'-> {bn[pick[t]]}')
