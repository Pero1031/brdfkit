"""
Surface and body of a coloured dielectric, split by polarization (Umov's law on one sample).

A pigmented plastic reflects in two ways: off its surface (Fresnel on tilted facets: white, and it keeps the polarization,
m11 = m00 in the table's frames) and from inside (the pigment's colour, almost fully depolarized). Where the pigment absorbs,
only the surface is left; where it does not, the body adds unpolarized light. So over wavelength, at one geometry,

    M(lambda) = M_surf + B(lambda) N_body + k c(lambda) e00

with M_surf constant, B the body's spectrum, N_body its normalized Mueller matrix (nearly e00 alone) and c the spectrum of the
instrument's unpolarized floor (docs/FLOOR.md). The polarized elements do not follow the colour, and the share of polarized
light, m11 / R with R = m00 - k c, falls as 1 / R: the darker the wavelength, the stronger the polarization (Umov 1905).

`split` works on the reductions of `brdfkit quirks` (Mueller matrices summed over phi_d in a theta_h range beside the lobe,
'mid', or far from it, 'far'), so it needs no pass over the table. The floor's spectrum c is m00 - m11 of a dark reference
(black glass) in the same region; k and a white remainder u are fitted on the bands where the body is black (`dark_max`):
m00 - m11 = k c + u there. Per theta_d it returns
    k, u        floor scale relative to the reference, and the unpolarized white remainder
    R           m00 - k c: the reflection with the floor removed [n_ch]
    surf        mean Mueller matrix over the dark bands: the surface term (its m00 still holds k c + u; use m11 as its level)
    body        slope of every element against R: the body's normalized Mueller matrix
    umov_slope  slope of log(m11 / R) against log R (-1 = Umov's law; NaN when R varies by less than 2x)
docs/UMOV.md.
"""
import numpy as np

REGIONS = ('mid', 'far')


def region_mean(red, region='mid'):
    """mean Mueller matrices of a quirks reduction in one theta_h region, [n_thd, n_ch, 4, 4]"""
    if region not in REGIONS:
        raise ValueError(f'region must be one of {REGIONS}')
    with np.errstate(invalid='ignore', divide='ignore'):
        return red[f'{region}_M'] / red[f'{region}_n'][:, None, None, None]


def fresnel_ratios(theta, n):
    """(m01 / m00, m22 / m00) of a Fresnel reflection off a dielectric of index n at the angle theta [rad], (s, p) frames with
    the table's signs (m01 > 0: s reflects more)"""
    ci = np.cos(theta); ct = np.sqrt(1 - (np.sin(theta) / n) ** 2)
    rs = (ci - n * ct) / (ci + n * ct); rp = (n * ci - ct) / (n * ci + ct)
    return (rs ** 2 - rp ** 2) / (rs ** 2 + rp ** 2), 2 * rs * rp / (rs ** 2 + rp ** 2)


def split(red, reference, region='mid', dark_max=540.0, theta_d_deg=range(25, 80, 5), lam_max=900.0, k_fixed=None):
    """-> (channels, list of dicts, one per theta_d). red, reference: reductions of `brdfkit quirks` (sample, dark reference)"""
    lam = np.asarray(red['channels'], float); td = np.degrees(np.asarray(red['theta_d'], float))
    if len(reference['channels']) != len(lam) or not np.allclose(reference['channels'], lam, atol=0.5) \
            or len(reference['theta_d']) != len(td):
        raise ValueError('the reference reduction is on another grid')
    M = region_mean(red, region); G = region_mean(reference, region)
    c = G[..., 0, 0] - G[..., 1, 1]                                                         # the floor's spectrum [T, C]
    dark = lam <= dark_max; use = lam <= lam_max
    if dark.sum() < 3:
        raise ValueError('fewer than 3 bands below dark_max')
    rows = []
    for t_deg in theta_d_deg:
        t = int(np.argmin(np.abs(td - t_deg))); m = M[t]
        if not np.isfinite(m).all() or not np.isfinite(c[t]).all():
            continue
        y = m[:, 0, 0] - m[:, 1, 1]
        if k_fixed is None:
            (k, u), *_ = np.linalg.lstsq(np.stack([c[t, dark], np.ones(dark.sum())], 1), y[dark], rcond=None)
        else:
            k, u = float(k_fixed), float(np.mean(y[dark] - k_fixed * c[t, dark]))
        R = m[:, 0, 0] - k * c[t]
        surf = m[dark].mean(0)
        B = R - R[dark].mean()
        body = np.array([[np.polyfit(B[use], m[use, i, j], 1)[0] for j in range(4)] for i in range(4)])
        pos = use & (R > 0) & (m[:, 1, 1] > 0)
        slope = np.nan
        if pos.sum() > 10 and R[pos].max() > 2 * R[pos].min():
            slope = float(np.polyfit(np.log(R[pos]), np.log(m[pos, 1, 1] / R[pos]), 1)[0])
        rows.append(dict(theta_d=float(td[t]), k=float(k), u=float(u), R=R, M=m, surf=surf, body=body, umov_slope=slope,
                         contrast=float(R[use].max() / R[use].min()) if R[use].min() > 0 else float('inf')))
    return lam, rows


def report(name, lam, rows, region='mid', n=None, say=print):
    say(f'== {name}: surface and body by polarization, theta_h {"beside the lobe" if region == "mid" else "far from the lobe"} ({region})')
    fr = ' (Fresnel)' if n else ''
    say(f'  theta_d  floor k  u/m11 | surface: m11 x1e3  m01/m11 m10/m11{fr}  m22/m11{fr}  m33/m11 | body: m11    m01    m10    m22    m33 | R max/min  Umov slope')
    for r in rows:
        s, b = r['surf'], r['body']
        f1 = f2 = ''
        if n:
            D, C = fresnel_ratios(np.radians(r['theta_d']), n); f1, f2 = f' ({D:5.3f})', f' ({C:+6.3f})'
        say(f'  {r["theta_d"]:5.1f}   {r["k"]:6.3f}  {r["u"] / s[1, 1]:6.2f} | {s[1, 1] * 1e3:14.3f}   {s[0, 1] / s[1, 1]:6.3f}  {s[1, 0] / s[1, 1]:6.3f}{f1}'
            f'   {s[2, 2] / s[1, 1]:+6.3f}{f2}   {s[3, 3] / s[1, 1]:+6.3f} |    {b[1, 1]:+6.3f} {b[0, 1]:+6.3f} {b[1, 0]:+6.3f} {b[2, 2]:+6.3f} {b[3, 3]:+6.3f} |'
            f' {r["contrast"]:8.1f}  {r["umov_slope"]:9.2f}')
    pick = [int(np.argmin(np.abs(lam - x))) for x in (438, 502, 550, 598, 646, 702, 798)]
    say('  degree of polarization under unpolarized light, m10 / R, per wavelength [nm]:')
    say('  theta_d ' + ' '.join(f'{lam[c]:7.0f}' for c in pick))
    for r in rows:
        with np.errstate(invalid='ignore', divide='ignore'):
            say(f'  {r["theta_d"]:5.1f}   ' + ' '.join(f'{r["M"][c, 1, 0] / r["R"][c]:7.3f}' for c in pick))
