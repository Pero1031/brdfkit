"""
Quirks of a measured table: one pass over the file, then a set of checks that name what a renderer or an analysis has to know.

`reduce` reads every phi_d slice once and keeps what the checks need (tens of MB for a 13 GB hpBRDF table):
  Y        [n_phi, n_thd, n_thh]       luminance of m00 (CIE Y weights over the bands)
  m00      [n_thd, n_thh, n_ch]        m00 summed over phi_d (with the count of measured bins)
  core_M   [n_thd, n_ch, 4, 4]         Mueller matrices summed over phi_d and theta_h <= core_deg (the specular core)
  mid_M    [n_thd, n_ch, 4, 4]         the same over theta_h in mid_deg (a baseline beside the lobe: the body of a diffuse material)
  far_M    [n_thd, n_ch, 4, 4]         the same over theta_h >= far_deg (far from the lobe)
`measure` turns them into the checks below, `report` prints them. Optics (a dielectric whose index is fitted from its own
specular polarization, or a metal's n, k) add the radiometric check and the polarization against Fresnel. docs/QUIRKS.md.

  inpainting    theta_d ranges where the table is constant along theta_d (filled in, not measured)
  smoothness    roughness of the specular value along theta_d (median |second difference| of its log)
  peak          flatness of the specular peak (value at theta_h 0.28 deg / value at 0) and the lobe's half width per theta_d
  far field     level, spectrum and polarization far from the lobe, and how its luminance depends on theta_d
  energy        lobe energy within theta_h < 10 deg above the local baseline, divided by Fresnel (1 = correctly scaled)
  spectrum      irregularity of the specular spectrum along theta_d, and whether the body's spectrum varies the same way
  phi_d         one-sided arms of the lobe (value at phi / value at -phi) and reciprocity (phi_d 0 against 180 at high theta_d)
  polarization  m01 (diattenuation) against m10 (polarizance) at the specular in (s, p) bases, and against Fresnel
  wiggle        bin-to-bin fluctuation of the lobe energy along theta_d after a smooth trend is removed
"""
import json
import warnings
import numpy as np
from ..table import BRDFTable
from . import factor as FA
from . import calibrate as CA
from .lu_chipman import Metal, film_rs_rp, fresnel_rs_rp, mueller_from_amplitudes

trapezoid = getattr(np, 'trapezoid', None) or np.trapz
KIND = 'brdfkit quirks reduction v1'


def reduce(table: BRDFTable, core_deg=0.55, mid_deg=(15.0, 25.0), far_deg=40.0, log=print):
    """one pass over the table -> dict of reduced arrays (see module doc)"""
    th = np.degrees(table.theta_h); lam = table.channels; yw = FA.luminance_weights(lam)
    nP, nD, nH, nC = table.n_phi, table.n_thd, table.n_thh, table.n_ch
    core = th <= core_deg + 1e-9; mid = (th >= mid_deg[0]) & (th <= mid_deg[1]); far = th >= far_deg
    Y = np.full((nP, nD, nH), np.nan, np.float32)
    m00 = np.zeros((nD, nH, nC)); m00_n = np.zeros((nD, nH))
    acc = {k: (np.zeros((nD, nC, 4, 4)), np.zeros(nD)) for k in ('core', 'mid', 'far')}
    masks = {'core': core, 'mid': mid, 'far': far}
    for i, sl in table.slices():
        ok = table.valid(sl)
        v = np.where(ok[..., None, None, None], sl, 0.0).astype(np.float64)
        Y[i] = np.where(ok, v[..., 0, 0] @ yw, np.nan)
        m00 += v[..., 0, 0]; m00_n += ok
        for k, m in masks.items():
            acc[k][0][...] += v[:, m].sum(1); acc[k][1][...] += ok[:, m].sum(1)
        if log and i % 60 == 0:
            log(f'  reduce {i}/{nP}')
    out = dict(kind=KIND, phi_d=table.phi_d, theta_d=table.theta_d, theta_h=table.theta_h, channels=lam, Y=Y, m00=m00, m00_n=m00_n,
               meta=json.dumps(dict(source=str(table.meta.get('source', '')), core_deg=core_deg, mid_deg=list(mid_deg), far_deg=far_deg)))
    for k, (s, n) in acc.items():
        out[f'{k}_M'] = s; out[f'{k}_n'] = n
    return out


def save(path, red):
    np.savez_compressed(path, **red)


def load(path):
    z = np.load(path, allow_pickle=False)
    d = {k: z[k] for k in z.files}
    d['meta'] = json.loads(str(d['meta']))
    return d


# ----------------------------------------------------------------------------- helpers
def _nanmean(x, axis):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(x, axis)


def _nanmedian(x, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmedian(x, axis)


def _at(arr_deg, value):
    return int(np.argmin(np.abs(arr_deg - value)))


def _norm(M, n):
    """summed Mueller matrices [T, C, 4, 4] / counts -> means, then divided by m00"""
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = M / n[:, None, None, None]
        return mean / mean[..., :1, :1]


def metal_fresnel(theta_d, channels, nk_path, film_nm=0.0, n_film=1.65):
    """(F [T, C], normalized Mueller [T, C, 4, 4]) of a metal (n, k table), optionally under a dielectric film"""
    met = Metal(nk_path); cos = np.cos(np.asarray(theta_d, float))
    F = np.zeros((len(cos), len(channels))); M = np.zeros((len(cos), len(channels), 4, 4))
    for c, lam in enumerate(channels):
        eta = met.eta(float(lam))
        rs, rp = film_rs_rp(cos, float(lam), eta, n_film, film_nm) if film_nm > 0 else fresnel_rs_rp(cos, eta)
        F[:, c] = CA.reflectance(rs, rp); M[:, c] = mueller_from_amplitudes(rs, rp)
    return F, M


# ----------------------------------------------------------------------------- the checks
def measure(red, optics=None, lobe_cut_deg=10.0, arm_theta_d=(25.0, 45.0), arm_theta_h=(6.0, 15.0)):
    """reduced arrays -> dict of checks. optics: None, ('dielectric',) or ('metal', nk_path, film_nm)"""
    with warnings.catch_warnings(), np.errstate(invalid='ignore', divide='ignore'):
        warnings.simplefilter('ignore', RuntimeWarning)   # theta_d bins without any measured value give NaN, on purpose
        return _measure(red, optics, lobe_cut_deg, arm_theta_d, arm_theta_h)


def _measure(red, optics, lobe_cut_deg, arm_theta_d, arm_theta_h):
    td =np.asarray(red['theta_d'], float); tdd = np.degrees(td); th = np.asarray(red['theta_h'], float); thd = np.degrees(th)
    ph = np.degrees(np.asarray(red['phi_d'], float)); lam = np.asarray(red['channels'], float); yw = FA.luminance_weights(lam)
    Y = np.asarray(red['Y'], np.float64)
    Ybar = _nanmean(Y, 0)                                                                # [T, H]
    c550 = _at(lam, 550.0); out = {}

    # inpainting: the specular value levels off toward both ends (a Gaussian fill blends into the measured range, so the edge is soft);
    # the theta_d below / above which it changes by less than 1e-4, 1e-3 and 1e-2 per bin
    v0 = Ybar[:, 0]
    with np.errstate(invalid='ignore', divide='ignore'):
        rel = np.abs(np.diff(v0)) / np.abs(v0[:-1])
    flat = {}
    for tol in (1e-4, 1e-3, 1e-2):
        moving = np.nonzero(rel > tol)[0]
        flat[tol] = (float(tdd[moving[0]]), float(tdd[moving[-1] + 1])) if len(moving) else (np.nan, np.nan)
    out['flat_ranges'] = flat
    lo, hi = flat[1e-3]
    in_rng = (tdd >= 20) & (tdd <= 79.5)                                                # the measured range (illumination 40-160 deg)
    out['measured_used'] = [float(tdd[in_rng][0]), float(tdd[in_rng][-1])]

    # smoothness of the specular value along theta_d
    lv = np.log(np.maximum(v0, 1e-12)); d2 = np.abs(lv[2:] - 2 * lv[1:-1] + lv[:-2])
    out['smooth_d2_measured'] = float(_nanmedian(d2[in_rng[1:-1]]))
    lowr = (tdd < lo - 1e-6) if np.isfinite(lo) else tdd < 15
    out['smooth_d2_inpainted'] = float(_nanmedian(d2[lowr[1:-1]])) if lowr[1:-1].any() else np.nan

    # peak flatness and lobe half width (above the baseline beside the lobe)
    mid = (thd >= 20) & (thd <= 30)
    base = _nanmedian(Ybar[:, mid], 1)
    out['peak_flat'] = float(_nanmedian((Ybar[:, _at(thd, 0.278)] / Ybar[:, 0])[in_rng]))
    hw = np.full(len(td), np.nan)
    for t in np.nonzero(in_rng)[0]:
        p = Ybar[t] - base[t]
        if not (p[0] > 0):
            continue
        k = np.nonzero(p < p[0] / 2)[0]
        if len(k) and k[0] > 0:
            a, b = k[0] - 1, k[0]
            hw[t] = thd[a] + (p[a] - p[0] / 2) / (p[a] - p[b]) * (thd[b] - thd[a])
    out['hwhm_deg'] = hw
    out['base_over_peak'] = float(_nanmedian((base / Ybar[:, 0])[in_rng]))

    # far field
    farN = _norm(red['far_M'], red['far_n'])
    far_m00 = red['far_M'][..., 0, 0] / red['far_n'][:, None]                          # [T, C]
    out['far_level_550'] = float(_nanmedian(far_m00[in_rng, c550]))
    with np.errstate(invalid='ignore', divide='ignore'):
        out['far_spec_414'] = float(_nanmedian((far_m00[:, 0] / far_m00[:, c550])[in_rng]))
        out['far_spec_950'] = float(_nanmedian((far_m00[:, -1] / far_m00[:, c550])[in_rng]))
    off = np.abs(farN[:, c550].reshape(len(td), 16)[:, 1:])
    out['far_pol_max_550'] = float(_nanmedian(np.nanmax(off, 1)[in_rng]))
    farY = far_m00 @ yw
    out['far_Y_vs_theta_d'] = farY / farY[_at(tdd, 30.0)]

    # lobe energy above the baseline, per theta_d (luminance), and the specular Mueller matrix above the body
    sel = thd <= lobe_cut_deg + 1e-9
    w = np.cos(th[sel]) * np.sin(th[sel])
    Q = trapezoid(np.nan_to_num(Ybar[:, sel] - base[:, None]) * w, th[sel], axis=1) * 2 * np.pi * 4 * np.cos(td) ** 2
    out['Q'] = Q
    coreM = red['core_M'] / red['core_n'][:, None, None, None]
    midM = red['mid_M'] / red['mid_n'][:, None, None, None]
    lobeM = coreM - midM                                                                # the body's Mueller matrix taken out
    with np.errstate(invalid='ignore', divide='ignore'):
        Msp = CA.table_to_sp(lobeM / lobeM[..., :1, :1], td)
    out['N01_550'] = Msp[:, c550, 0, 1]; out['N10_550'] = Msp[:, c550, 1, 0]; out['N02_550'] = Msp[:, c550, 0, 2]

    F = None
    if optics and optics[0] == 'dielectric':
        rows = in_rng & (tdd >= 25) & (tdd <= 75)
        n, n_rms = CA.fit_n(Msp, td, rows)
        F, Mf = CA.dielectric(td, n)
        out['n_fit_median'], out['n_fit_rms'] = float(np.median(n)), float(np.median(n_rms))
    elif optics and optics[0] == 'metal':
        F, Mf = metal_fresnel(td, lam, optics[1], optics[2] if len(optics) > 2 else 0.0)
    if F is not None:
        FY = F @ yw
        out['energy_over_fresnel'] = np.where(in_rng, Q / FY, np.nan)
        out['N01_fresnel_550'] = Mf[:, c550, 0, 1]

    # spectrum: the specular core's spectrum (body removed) and the far field's, both normalized per theta_d
    spec = np.maximum(lobeM[..., 0, 0], 1e-12); spec = spec / (spec @ yw)[:, None]
    farS = np.maximum(far_m00, 1e-12); farS = farS / (farS @ yw)[:, None]
    def irregular(S):
        L = np.log(S[in_rng]); r = L[2:] - 2 * L[1:-1] + L[:-2]
        return float(np.median(np.sqrt(np.mean(r ** 2, 0))))
    out['spec_irregularity'] = irregular(spec); out['far_spec_irregularity'] = irregular(farS)
    rs = np.log(spec[in_rng]) - np.log(spec[in_rng]).mean(0); rf = np.log(farS[in_rng]) - np.log(farS[in_rng]).mean(0)
    out['spec_far_corr'] = float(np.corrcoef(rs.ravel(), rf.ravel())[0, 1])

    # phi_d arms and reciprocity
    hsel = (thd >= arm_theta_h[0]) & (thd <= arm_theta_h[1])
    mirror = np.array([_at(ph, -p) for p in ph])
    arms = {}
    for t0 in arm_theta_d:
        A = _nanmean(Y[:, _at(tdd, t0)][:, hsel], 1)
        ratio = A / A[mirror]
        k = int(np.nanargmax(ratio))
        arms[float(t0)] = dict(phi=float(ph[k]), arm_over_mirror=float(ratio[k]), arm_over_median=float(np.nanmax(A) / _nanmedian(A)))
    out['arms'] = arms
    p0, p180 = _at(ph, 0.0), _at(ph, 180.0)
    hr = (thd >= 1) & (thd <= 10); trs = (tdd >= 65) & in_rng
    with np.errstate(invalid='ignore', divide='ignore'):
        out['reciprocity_0_over_180'] = float(_nanmedian(Y[p0][trs][:, hr] / Y[p180][trs][:, hr]))
        lowd = (tdd >= 25) & (tdd <= 45)
        out['reciprocity_0_over_180_low'] = float(_nanmedian(Y[p0][lowd][:, hr] / Y[p180][lowd][:, hr]))

    # wiggle of the lobe energy along theta_d
    x = tdd[in_rng]; lq = np.log(np.maximum(Q[in_rng], 1e-12))
    out['wiggle_rms'] = float(np.std(lq - np.polyval(np.polyfit(x, lq, 3), x)))
    out['theta_d_deg'] = tdd
    return out


def report(name, m, say=print):
    tdd = m['theta_d_deg']; pick = [t for t in (20, 30, 40, 50, 60, 70, 75) if m['measured_used'][0] <= t <= m['measured_used'][1]]
    idx = [int(np.argmin(np.abs(tdd - t))) for t in pick]
    f = lambda a: '  '.join(f'{a[i]:6.3f}' if np.isfinite(a[i]) else '     -' for i in idx)   # noqa: E731
    say(f'== {name}')
    say('  inpainting: the specular value changes by less than ' + ', '.join(
        f'{tol:.0e} per bin below {a:.0f} / above {b:.0f} deg' for tol, (a, b) in m['flat_ranges'].items())
        + f'; checks use {m["measured_used"][0]:.0f}-{m["measured_used"][1]:.0f} deg')
    say(f'  smoothness: median |2nd difference of log specular| {m["smooth_d2_measured"]:.4f} (measured), {m["smooth_d2_inpainted"]:.1e} (inpainted)')
    say(f'  peak: value at 0.28 deg / at 0 = {m["peak_flat"]:.3f}; baseline beside the lobe / peak = {m["base_over_peak"]:.2e}')
    say('  theta_d               ' + '  '.join(f'{t:6d}' for t in pick))
    say('  lobe half width [deg] ' + f(m['hwhm_deg']))
    if 'energy_over_fresnel' in m:
        say('  lobe energy / Fresnel ' + f(m['energy_over_fresnel']))
    say('  far luminance / 30deg ' + f(m['far_Y_vs_theta_d']))
    say('  N01 (s-p) at 550 nm   ' + f(m['N01_550']))
    say('  N10 (column)          ' + f(m['N10_550']))
    if 'N01_fresnel_550' in m:
        say('  N01 Fresnel           ' + f(m['N01_fresnel_550']))
    if 'n_fit_median' in m:
        say(f'  index fitted from the specular polarization: n = {m["n_fit_median"]:.3f} (residual {m["n_fit_rms"]:.3f})')
    say(f'  far field: level {m["far_level_550"]:.4f} at 550 nm, spectrum 414 / 550 = {m["far_spec_414"]:.2f}, 950 / 550 = {m["far_spec_950"]:.2f},'
        f' max |m_ij / m00| {m["far_pol_max_550"]:.3f}')
    say(f'  spectrum: irregularity along theta_d (RMS 2nd difference of log) specular {m["spec_irregularity"]:.3f}, far field'
        f' {m["far_spec_irregularity"]:.3f}; correlation of their deviations {m["spec_far_corr"]:+.2f}')
    for t0, a in m['arms'].items():
        say(f'  phi_d arm at theta_d {t0:.0f}: phi {a["phi"]:+.0f} deg, {a["arm_over_mirror"]:.2f} x its mirror, {a["arm_over_median"]:.2f} x the median')
    say(f'  reciprocity: phi_d 0 / 180 = {m["reciprocity_0_over_180"]:.2f} at theta_d >= 65, {m["reciprocity_0_over_180_low"]:.2f} at 25-45 (theta_h 1-10)')
    say(f'  wiggle of the lobe energy along theta_d: {m["wiggle_rms"]:.3f} (RMS of log after a cubic trend)')
