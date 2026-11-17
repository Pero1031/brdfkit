"""
Check an hpBRDF table against a reference with known reflection, black glass, and correct the radiometric error it reveals.

A polished dielectric reflects as Fresnel with one real index n(lambda). From the factor model of the glass
(`brdfkit factor`, floor removed; m00 = A x S, m_ij/m00 = N, docs/FACTOR.md):

  1. Polarization. N(theta_d, lambda) is rotated into the (s, p) bases at the specular geometry and fitted with the Fresnel
     Mueller matrix of a dielectric (m01, m10, m22, m33 over the measured theta_d range) -> n(lambda). The residual and the
     Lu-Chipman diattenuation and retardance against the fit show the polarimetric error of the instrument at the specular.
  2. Reflectance. The energy of the lobe around the specular direction,
         Q(theta_d) = 2 pi 4 cos^2(theta_d) int_0^cutoff  Abar(theta_d, theta_h) cos(theta_h) sin(theta_h) d theta_h
     (Abar: A averaged over phi_d), is the directional albedo of a narrow lobe of a correctly scaled BRDF (theta_i = theta_o
     = theta_d, D normalized), so for the glass it should equal its Fresnel reflectance F_Y(theta_d). K = Q / F_Y is the
     radiometric error of the table at that theta_d; in the hpBRDF instrument one theta_d is one illumination angle
     (theta_d = half the light-camera angle). K is smoothed along theta_d and copied from the edges outside the measured range.
  3. Spectral shape. S(theta_d, lambda) against F(theta_d, lambda) / F_Y (1 = the shape of Fresnel; glass is nearly flat).

`apply(model, calib)` divides A of any factor model by K(theta_d): luminance only, S and N are kept (the spectral and
polarimetric errors of the glass do not repeat on other materials, docs/CALIBRATION.md). `check(model, calib, nk)` compares
the directional albedo of a metal with its Fresnel reflectance before and after. Calibration file (.npz): kind =
'brdfkit-calibration', theta_d, K, K_raw, Q, F_Y, n, n_rms, channels, spectral_error, D_meas, D_fresnel, R_meas, R_fresnel,
m02, cutoff_deg, smooth_deg, measured_theta_d, source.
"""
import json
import warnings
import numpy as np
from scipy.ndimage import gaussian_filter1d
from .. import geometry as G
from . import factor as FA
from .lu_chipman import Metal, film_rs_rp, fresnel_rs_rp, lu_chipman, mueller_from_amplitudes, specular_directions

KIND = 'brdfkit-calibration'
FIT_ELEMENTS = ((0, 1), (1, 0), (2, 2), (3, 3))


# ----------------------------------------------------------------------------- frames
def _frames(theta_deg):
    tl, tv = specular_directions(theta_deg)
    g = G.geometry(tl, tv)
    return tl, tv, g['x_in'], g['x_out'], G.s_axis(tl, tv)


def table_to_sp(M, theta_d):
    """normalized table-frame Mueller matrices [n_thd, n_ch, 4, 4] at the specular geometry -> (s, p) bases (NaN at 0, 90 deg)"""
    th = np.degrees(np.asarray(theta_d, float)); out = np.full(M.shape, np.nan)
    ok = (th > 0.5) & (th < 89.5)
    tl, tv, xi, xo, s = _frames(th[ok])
    for c in range(M.shape[1]):
        out[ok, c] = G.rotate_mueller_basis(M[ok, c], -tl, xi, s, tv, xo, s)
    return out


def sp_to_table(M, theta_d):
    """inverse of `table_to_sp`"""
    th = np.degrees(np.asarray(theta_d, float)); out = np.full(M.shape, np.nan)
    ok = (th > 0.5) & (th < 89.5)
    tl, tv, xi, xo, s = _frames(th[ok])
    for c in range(M.shape[1]):
        out[ok, c] = G.rotate_mueller_basis(M[ok, c], -tl, s, xi, tv, s, xo)
    return out


def specular_sp(model):
    """the polarization N of a factor model as full normalized Mueller matrices in (s, p) bases, [n_thd, n_ch, 4, 4]"""
    N = np.asarray(model['N'], float)
    M = np.concatenate([np.ones(N.shape[:2] + (1,)), N], -1).reshape(N.shape[:2] + (4, 4))
    return table_to_sp(M, model['theta_d'])


# ----------------------------------------------------------------------------- Fresnel references
def reflectance(rs, rp):
    return 0.5 * (np.abs(rs) ** 2 + np.abs(rp) ** 2)


def dielectric(theta_d, n):
    """(F [T, C], normalized Mueller [T, C, 4, 4]) of a dielectric with index n [C] at angles theta_d [T]"""
    cos = np.cos(np.asarray(theta_d, float))[:, None]
    rs, rp = fresnel_rs_rp(np.broadcast_to(cos, (len(cos), len(n))), np.asarray(n, float)[None, :] + 0j)
    return reflectance(rs, rp), mueller_from_amplitudes(rs, rp)


def metal_reflectance(theta_d, channels, nk_path, film_nm=0.0, n_film=1.65):
    """F [T, C] of a metal (refractiveindex.info n, k table), optionally under one dielectric film"""
    met = Metal(nk_path); cos = np.cos(np.asarray(theta_d, float))
    F = np.zeros((len(cos), len(channels)))
    for c, lam in enumerate(channels):
        eta = met.eta(float(lam))
        F[:, c] = reflectance(*(film_rs_rp(cos, float(lam), eta, n_film, film_nm) if film_nm > 0 else fresnel_rs_rp(cos, eta)))
    return F


# ----------------------------------------------------------------------------- estimation
def lobe_energy(model, cutoff_deg=10.0):
    """Q(theta_d) [n_thd]: the directional albedo of the lobe near the specular (luminance, since S has Y = 1)"""
    A = np.asarray(model['A'], np.float64); td = np.asarray(model['theta_d'], float); th = np.asarray(model['theta_h'], float)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        Abar = np.nanmean(A, 0)
    sel = np.degrees(th) <= cutoff_deg + 1e-9
    w = np.cos(th[sel]) * np.sin(th[sel])
    return np.trapz(np.nan_to_num(Abar[:, sel]) * w, th[sel], axis=1) * 2 * np.pi * 4 * np.cos(td) ** 2


def fit_n(Msp, theta_d, rows, n_range=(1.3, 2.2)):
    """per channel: the real index whose Fresnel Mueller matrix best matches m01, m10, m22, m33 over `rows` -> (n [C], rms [C])"""
    T = np.asarray(theta_d, float)[rows]; Y = Msp[rows]
    def cost(grid, c):
        _, Mm = dielectric(T, grid)                                                   # [T, G, 4, 4]
        e = sum((Mm[:, :, i, j] - Y[:, c, i, j][:, None]) ** 2 for i, j in FIT_ELEMENTS)
        return np.nanmean(e, 0) / len(FIT_ELEMENTS)
    n = np.zeros(Msp.shape[1]); rms = np.zeros(Msp.shape[1])
    coarse = np.arange(n_range[0], n_range[1], 0.01)
    for c in range(Msp.shape[1]):
        g0 = coarse[np.argmin(cost(coarse, c))]
        fine = np.arange(g0 - 0.01, g0 + 0.01, 0.0002); e = cost(fine, c)
        n[c] = fine[np.argmin(e)]; rms[c] = np.sqrt(e.min())
    return n, rms


def estimate(glass, cutoff_deg=10.0, smooth_deg=3.0, measured=None):
    """calibration dict from the factor model of a polished dielectric (black glass)"""
    td = np.asarray(glass['theta_d'], float); tdd = np.degrees(td); lam = np.asarray(glass['channels'], float)
    measured = tuple(np.asarray(glass['measured_theta_d'], float)) if measured is None else tuple(measured)
    meas = (tdd >= measured[0] - 1e-6) & (tdd <= measured[1] + 1e-6)
    Msp = specular_sp(glass)
    n, n_rms = fit_n(Msp, td, meas)
    F, Mf = dielectric(td, n)
    yw = FA.luminance_weights(lam); FY = F @ yw
    Q = lobe_energy(glass, cutoff_deg)
    K_raw = np.where(meas, Q / FY, np.nan)
    K = K_raw.copy()
    if smooth_deg > 0:
        idx = np.nonzero(meas)[0]
        K[idx] = gaussian_filter1d(K_raw[idx], smooth_deg / np.median(np.diff(tdd)), mode='nearest')
    K = FA._extend_edges(K, meas)
    S = np.asarray(glass['S'], float)
    E = S / (F / FY[:, None])
    D_meas = np.full(F.shape, np.nan); R_meas = np.full(F.shape, np.nan); D_f = np.full(F.shape, np.nan); R_f = np.full(F.shape, np.nan)
    for t in np.nonzero(meas)[0]:
        for c in range(len(lam)):
            d = lu_chipman(Msp[t, c]); f = lu_chipman(Mf[t, c])
            if d is not None:
                D_meas[t, c], R_meas[t, c] = d['D'], d['R']
            if f is not None:
                D_f[t, c], R_f[t, c] = f['D'], f['R']
    return dict(kind=KIND, theta_d=td, K=K, K_raw=K_raw, Q=Q, F_Y=FY, n=n, n_rms=n_rms, channels=lam, spectral_error=E,
                D_meas=D_meas, D_fresnel=D_f, R_meas=R_meas, R_fresnel=R_f, m02=Msp[..., 0, 2], cutoff_deg=float(cutoff_deg),
                smooth_deg=float(smooth_deg), measured_theta_d=np.array(measured, float),
                source=str(glass.get('meta', {}).get('source', '')))


def save(path, calib):
    np.savez(path, **calib)
    return path


def load(path):
    with np.load(path, allow_pickle=False) as z:
        c = {k: z[k] for k in z.files}
    if str(c.get('kind')) != KIND:
        raise ValueError(f'{path}: not a brdfkit calibration file')
    for k in ('cutoff_deg', 'smooth_deg'):
        c[k] = float(c[k])
    c['source'] = str(c['source'])
    return c


def K_at(calib, theta_d):
    """K interpolated (in degrees) at angles theta_d [rad]"""
    return np.interp(np.degrees(np.asarray(theta_d, float)), np.degrees(calib['theta_d']), calib['K'])


def apply(model, calib, keep_level_deg=None):
    """factor model with A divided by K(theta_d). keep_level_deg: rescale so that the level at that theta_d is unchanged
    (corrects only the theta_d dependence, not the absolute scale)"""
    K = K_at(calib, model['theta_d'])
    scale = 1.0 / K
    if keep_level_deg is not None:
        scale = scale * float(K_at(calib, np.radians([keep_level_deg]))[0])
    out = dict(model)
    out['A'] = (np.asarray(model['A'], np.float32) * scale[None, :, None].astype(np.float32))
    out['meta'] = dict(model.get('meta', {}), radiometric_calibration=dict(source=calib['source'], cutoff_deg=calib['cutoff_deg'],
                       smooth_deg=calib['smooth_deg'], keep_level_deg=keep_level_deg))
    return out


def check(model, calib, nk_path, film_nm=0.0, n_film=1.65, angles=(25, 35, 45, 55, 65, 75), n_mc=1_000_000, keep_level_deg=None):
    """directional albedo / Fresnel reflectance of a metal, before and after `apply` -> list of (angle, albedo, F, ratio, ratio after)"""
    from .latent_edit import albedo
    lam = np.asarray(model['channels'], float); yw = FA.luminance_weights(lam)
    after = apply(model, calib, keep_level_deg)
    rows = []
    for a in angles:
        F = float(metal_reflectance(np.radians([a]), lam, nk_path, film_nm, n_film)[0] @ yw)
        b, c = albedo(model, a, n=n_mc), albedo(after, a, n=n_mc)
        rows.append((float(a), b, F, b / F, c / F))
    return rows


# ----------------------------------------------------------------------------- report
def report(calib, check_rows=None, check_name='', say=print):
    td = np.degrees(calib['theta_d']); lam = calib['channels']; meas = np.isfinite(calib['K_raw'])
    say(f'calibration from {calib["source"]}  (lobe energy within theta_h < {calib["cutoff_deg"]:g} deg, K smoothed {calib["smooth_deg"]:g} deg, '
        f'measured theta_d {calib["measured_theta_d"][0]:g}-{calib["measured_theta_d"][1]:g})')
    n, r = calib['n'], calib['n_rms']
    say(f'index of the glass from its polarization: n = {np.median(n):.3f} (range {n.min():.3f}-{n.max():.3f}); '
        f'fit residual (RMS of m01, m10, m22, m33) median {np.median(r):.3f}, max {r.max():.3f}')
    say('  ' + '  '.join(f'{lam[c]:.0f} nm {n[c]:.3f}' for c in range(0, len(lam), max(1, len(lam) // 8))))
    say('radiometric error K = lobe energy / Fresnel reflectance (1 = correct; the same K on two materials = the instrument):')
    say('  theta_d  ' + ' '.join(f'{t:6.0f}' for t in range(20, 81, 5)))
    rows = [int(np.argmin(abs(td - t))) for t in range(20, 81, 5)]
    say('  K raw    ' + ' '.join(f'{calib["K_raw"][k]:6.3f}' for k in rows))
    say('  K used   ' + ' '.join(f'{calib["K"][k]:6.3f}' for k in rows))
    say('  cos^2    ' + ' '.join(f'{np.cos(calib["theta_d"][k]) ** 2:6.3f}' for k in rows))
    cs = [int(np.argmin(abs(lam - l))) for l in (414, 462, 550, 702, 798, 942)]
    say('spectral shape error of the glass S / (F / F_Y) (1 = Fresnel shape):')
    say('  theta_d  ' + ' '.join(f'{lam[c]:6.0f}' for c in cs))
    for t in (20, 30, 40, 50, 60, 70, 80):
        k = int(np.argmin(abs(td - t)))
        say(f'  {t:5d}    ' + ' '.join(f'{calib["spectral_error"][k, c]:6.2f}' for c in cs))
    c550 = int(np.argmin(abs(lam - 550)))
    say(f'polarization of the glass at {lam[c550]:.0f} nm against Fresnel(n): diattenuation D, retardance R [deg], m02 (ideal 0)')
    for t in (20, 25, 30, 40, 50, 60, 70, 80):
        k = int(np.argmin(abs(td - t)))
        if meas[k]:
            say(f'  {t:5d}   D {calib["D_meas"][k, c550]:.3f} / {calib["D_fresnel"][k, c550]:.3f}   '
                f'R {np.degrees(calib["R_meas"][k, c550]):6.1f} / {np.degrees(calib["R_fresnel"][k, c550]):6.1f}   m02 {calib["m02"][k, c550]:+.3f}')
    if check_rows:
        say(f'check on {check_name}: directional albedo (Monte Carlo) / Fresnel reflectance, before -> after the correction')
        for a, b, F, r0, r1 in check_rows:
            say(f'  theta_i {a:4.0f}   albedo {b:.3f}   Fresnel {F:.3f}   ratio {r0:.3f} -> {r1:.3f}')


def plot(calib, path, check_rows=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    td = np.degrees(calib['theta_d']); lam = calib['channels']; meas = np.isfinite(calib['K_raw'])
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    a = ax[0, 0]
    a.plot(td[meas], calib['K_raw'][meas], '.', color='grey', label='K raw (black glass)')
    a.plot(td, calib['K'], '-', color='C0', label='K used')
    a.plot(td, 0.55 * np.cos(np.radians(td)) ** 2, ':', color='k', label='0.55 cos²θ_d (reference)')
    if check_rows:
        a.plot([r[0] for r in check_rows], [r[3] for r in check_rows], 's', color='C3', label='metal: albedo / Fresnel (check)')
    a.set_xlabel('θ_d [deg]'); a.set_ylabel('lobe energy / Fresnel'); a.set_title('radiometric error of the table'); a.legend(fontsize=8)
    a = ax[0, 1]
    a.plot(lam, calib['n'], '.-'); a.set_xlabel('wavelength [nm]'); a.set_ylabel('n'); a.set_title('index of the glass from its polarization')
    c550 = int(np.argmin(abs(lam - 550)))
    a = ax[1, 0]
    a.plot(td[meas], calib['D_meas'][meas, c550], '.', label='measured'); a.plot(td[meas], calib['D_fresnel'][meas, c550], '-', label='Fresnel(n)')
    a.plot(td[meas], calib['m02'][meas, c550], 'x', ms=4, label='m02 (ideal 0)')
    a.set_xlabel('θ_d [deg]'); a.set_title(f'glass diattenuation at {lam[c550]:.0f} nm'); a.legend(fontsize=8)
    a = ax[1, 1]
    im = a.imshow(np.log2(calib['spectral_error'][meas]), aspect='auto', origin='lower', cmap='RdBu_r', vmin=-1, vmax=1,
                  extent=[lam[0], lam[-1], td[meas][0], td[meas][-1]])
    fig.colorbar(im, ax=a, label='log2 S / Fresnel shape'); a.set_xlabel('wavelength [nm]'); a.set_ylabel('θ_d [deg]')
    a.set_title('spectral shape error of the glass')
    fig.tight_layout(); fig.savefig(path, dpi=100); plt.close(fig)
    return path
