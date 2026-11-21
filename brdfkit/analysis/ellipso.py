"""
The specular reflection of a table as an ellipsometer: which material is it?

A reflection off a smooth isotropic surface is described by rho = rp / rs = tan(Psi) exp(i Delta). In the normalized Mueller
matrix of the specular geometry, in (s, p) bases, the phase Delta sits in m22, m33, m23, m32 alone:

    (m22 + m33) / 2 = sin(2 Psi) cos(Delta),   (m23 - m32) / 2 = sin(2 Psi) sin(Delta)      (signs as lu_chipman.mueller_from_amplitudes)

so Delta = atan2 of the two does not need the first row and column, which the hpBRDF tables get wrong on metals
(docs/QUIRKS.md, 9), nor the absolute scale, nor is it changed by depolarization, which only shrinks both terms. Delta is
measured at every theta_d (the angle of incidence at the specular) and every wavelength, and it falls from 180 deg at normal
incidence through 90 deg at the principal angle: a curve that depends on n and k.

`delta(reduction)` reads it from the specular core of a `brdfkit quirks` reduction. `fit_nk` fits n, k per wavelength to
Delta over theta_d (pseudo optical constants: those of a bare surface with that phase). `match` scores candidate materials
(n, k tables), bare or under one transparent film whose thickness is fitted, by the RMS difference of Delta in degrees.
docs/ELLIPSO.md.
"""
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from . import calibrate as CA
from .lu_chipman import Metal, film_rs_rp, fresnel_rs_rp

THETA_D = (25.0, 75.0)        # theta_d range used [deg]: inside the measured one, where the lobe still resolves


def delta_of_mueller(Msp):
    """(Delta [deg] in (-180, 180], q = sin(2 Psi) x polarization kept) of normalized Mueller matrices in (s, p) bases"""
    s = 0.5 * (Msp[..., 2, 3] - Msp[..., 3, 2]); c = 0.5 * (Msp[..., 2, 2] + Msp[..., 3, 3])
    return np.degrees(np.arctan2(s, c)), np.hypot(s, c)


def delta(red, theta_d=THETA_D):
    """specular core of a quirks reduction -> dict(theta_d [deg], channels, delta [T, C], q [T, C]) over the theta_d range"""
    td = np.asarray(red['theta_d'], float); tdd = np.degrees(td)
    with np.errstate(invalid='ignore', divide='ignore'):
        M = red['core_M'] / red['core_n'][:, None, None, None]
        Msp = CA.table_to_sp(M / M[..., :1, :1], td)
    rows = (tdd >= theta_d[0] - 1e-6) & (tdd <= theta_d[1] + 1e-6)
    d, q = delta_of_mueller(Msp[rows])
    return dict(theta_d=tdd[rows], channels=np.asarray(red['channels'], float), delta=d, q=q)


def model_delta(theta_d_deg, eta, lam=None, film_nm=0.0, n_film=1.5):
    """Delta [deg] of a surface of complex index eta (scalar or [C]) at angles theta_d [T] -> [T] or [T, C].
    film_nm > 0: under one transparent film (needs lam [nm], scalar or [C])"""
    cos = np.cos(np.radians(np.asarray(theta_d_deg, float)))
    eta = np.asarray(eta, complex)
    if eta.ndim == 0:
        rs, rp = film_rs_rp(cos, float(lam), eta, n_film, film_nm) if film_nm > 0 else fresnel_rs_rp(cos, eta)
        return np.degrees(np.angle(rs * np.conj(rp)))
    lam = np.broadcast_to(np.asarray(lam, float), eta.shape) if lam is not None else [None] * len(eta)
    return np.stack([model_delta(theta_d_deg, e, l, film_nm, n_film) for e, l in zip(eta, lam)], 1)


def _wrap(d):
    return (d + 180.0) % 360.0 - 180.0


def fit_nk(theta_d_deg, delta_deg):
    """pseudo optical constants: per wavelength the n, k of a bare surface whose Delta(theta_d) fits best.
    delta_deg [T, C] -> (n [C], k [C], rms [C] in degrees)"""
    T, C = delta_deg.shape
    n = np.zeros(C); k = np.zeros(C); rms = np.zeros(C)
    grid = [(a, b) for a in (0.2, 0.7, 1.5, 2.5, 3.5) for b in (0.5, 1.5, 3.0, 5.0, 7.0)]
    for c in range(C):
        y = delta_deg[:, c]; ok = np.isfinite(y)
        res = lambda p: _wrap(model_delta(theta_d_deg[ok], complex(abs(p[0]), abs(p[1]))) - y[ok])     # noqa: E731
        p0 = min(grid, key=lambda p: np.sum(res(p) ** 2))
        s = least_squares(res, p0)
        n[c], k[c] = abs(s.x[0]), abs(s.x[1]); rms[c] = np.sqrt(np.mean(s.fun ** 2))
    return n, k, rms


def match(meas, candidates, lam_range=(430.0, 900.0), film_max=200.0, n_film=1.5):
    """meas: dict of `delta`; candidates: {name: path of an n, k table (.yml)} -> list of dicts sorted by the best RMS:
    name, rms_bare [deg], film_nm, rms_film [deg] (one film thickness for all wavelengths, fitted in 0..film_max)"""
    lam = meas['channels']; td = meas['theta_d']; out = []
    for name, path in candidates.items():
        met = Metal(path)
        use = (lam >= max(lam_range[0], met.wl.min())) & (lam <= min(lam_range[1], met.wl.max()))
        if use.sum() < 3:
            continue
        eta = np.array([met.eta(l) for l in lam[use]]); y = meas['delta'][:, use]

        def rms(d):
            return float(np.sqrt(np.nanmean(_wrap(model_delta(td, eta, lam[use], d, n_film) - y) ** 2)))
        coarse = np.arange(0.0, film_max + 1e-9, 5.0); r = [rms(d) for d in coarse]; j = int(np.argmin(r))
        best = minimize_scalar(rms, bounds=(max(coarse[j] - 5, 0), min(coarse[j] + 5, film_max)), method='bounded')
        out.append(dict(name=name, rms_bare=r[0], film_nm=float(best.x), rms_film=float(best.fun), n_bands=int(use.sum()),
                        lam_range=(float(lam[use][0]), float(lam[use][-1]))))
    return sorted(out, key=lambda x: x['rms_film'])


def principal_angle(theta_d_deg, delta_deg):
    """theta_d [deg] at which |Delta| passes 90 deg, per wavelength (NaN if it does not inside the range)"""
    a = np.abs(delta_deg); out = np.full(a.shape[1], np.nan)
    for c in range(a.shape[1]):
        j = np.nonzero((a[:-1, c] >= 90) & (a[1:, c] < 90))[0]
        if len(j):
            i = j[0]; out[c] = theta_d_deg[i] + (a[i, c] - 90) / (a[i, c] - a[i + 1, c]) * (theta_d_deg[i + 1] - theta_d_deg[i])
    return out


def report(name, meas, ranking, nk=None, say=print):
    lam = meas['channels']; td = meas['theta_d']
    pick = [int(np.argmin(np.abs(lam - x))) for x in (462, 502, 550, 598, 646, 702, 798, 902)]
    pa = principal_angle(td, meas['delta'])
    say(f'== {name}: phase of the specular reflection, theta_d {td[0]:.0f}-{td[-1]:.0f} deg')
    say('  wavelength [nm]        ' + ' '.join(f'{lam[c]:6.0f}' for c in pick))
    say('  principal angle [deg]  ' + ' '.join(f'{pa[c]:6.1f}' for c in pick))
    say('  q at 45 deg            ' + ' '.join(f'{meas["q"][int(np.argmin(np.abs(td - 45))), c]:6.2f}' for c in pick))
    if nk is not None:
        n, k, rms = nk
        say('  pseudo n               ' + ' '.join(f'{n[c]:6.2f}' for c in pick))
        say('  pseudo k               ' + ' '.join(f'{k[c]:6.2f}' for c in pick))
        say('  fit RMS [deg]          ' + ' '.join(f'{rms[c]:6.1f}' for c in pick))
    say('  candidate                      RMS of Delta bare [deg]   best film [nm]   RMS with it [deg]   range [nm]')
    for r in ranking:
        say(f'  {r["name"]:30s} {r["rms_bare"]:10.1f} {r["film_nm"]:20.0f} {r["rms_film"]:16.1f}          {r["lam_range"][0]:.0f}-{r["lam_range"][1]:.0f}')


# ----------------------------------------------------------------------------- phase + colour + level
TOLERANCE = dict(delta=3.0, shape=0.15, level=0.15)      # degrees; log units of the spectral shape; log units of the reflectance


def specular_reflectance(model, calib, theta_d=(25.0, 50.0)):
    """from a factor model and a calibration built with the same settings: dict(theta_d [deg], level [T] the luminance
    reflectance of the lobe (lobe energy / K), S [T, C] its spectrum normalized to luminance 1, channels)"""
    from .clean import check_calibration
    check_calibration(model.get('meta', {}), calib)
    tdd = np.degrees(np.asarray(model['theta_d'], float)); rows = (tdd >= theta_d[0] - 1e-6) & (tdd <= theta_d[1] + 1e-6)
    done = bool(model.get('meta', {}).get('radiometric_calibration'))                       # already divided by K (edit --calib)
    level = (CA.lobe_energy(model, calib['cutoff_deg']) / (1.0 if done else CA.K_at(calib, model['theta_d'])))[rows]
    return dict(theta_d=tdd[rows], level=level, S=np.asarray(model['S'], float)[rows], channels=np.asarray(model['channels'], float))


def _model(theta_d_deg, lam, met, film_nm, n_film):
    cos = np.cos(np.radians(theta_d_deg)); R = np.zeros((len(theta_d_deg), len(lam))); D = np.zeros_like(R)
    for c, l in enumerate(lam):
        rs, rp = film_rs_rp(cos, float(l), met.eta(l), n_film, film_nm) if film_nm > 0 else fresnel_rs_rp(cos, met.eta(l))
        R[:, c] = 0.5 * (abs(rs) ** 2 + abs(rp) ** 2); D[:, c] = np.degrees(np.angle(rs * np.conj(rp)))
    return R, D


def misfit(meas, refl, met, film_nm=0.0, n_film=1.5, lam_range=(430.0, 800.0), tol=TOLERANCE):
    """(combined, Delta RMS [deg], shape RMS [log], level RMS [log], model luminance reflectance) of one candidate.
    combined: RMS of the three divided by their tolerances (about 1 = as good as the data allow)"""
    from .factor import luminance_weights
    lam = meas['channels']; use = (lam >= lam_range[0]) & (lam <= lam_range[1]); yw = luminance_weights(lam)
    _, D = _model(meas['theta_d'], lam[use], met, film_nm, n_film)
    e_d = float(np.sqrt(np.nanmean(_wrap(D - meas['delta'][:, use]) ** 2)))
    R, _ = _model(refl['theta_d'], lam, met, film_nm, n_film); Y = R @ yw
    e_s = float(np.sqrt(np.mean(np.log((R / Y[:, None])[:, use] / refl['S'][:, use]) ** 2)))
    e_l = float(np.sqrt(np.mean(np.log(Y / refl['level']) ** 2)))
    return float(np.sqrt(((e_d / tol['delta']) ** 2 + (e_s / tol['shape']) ** 2 + (e_l / tol['level']) ** 2) / 3)), e_d, e_s, e_l, float(Y.mean())


def identify(meas, refl, candidates, film_max=120.0, n_film=1.5, lam_range=(430.0, 800.0)):
    """rank candidates {name: n, k table path} by `misfit`, each under its best film thickness in 0..film_max.
    -> list of dicts (name, film_nm, combined, delta, shape, level, reflectance, bare: the same without a film), best first"""
    out = []
    for name, path in candidates.items():
        met = Metal(path)
        f = lambda d: misfit(meas, refl, met, d, n_film, lam_range)[0]                 # noqa: E731
        coarse = np.arange(0.0, film_max + 1e-9, 5.0); r = [f(d) for d in coarse]; j = int(np.argmin(r))
        b = minimize_scalar(f, bounds=(max(coarse[j] - 5, 0), min(coarse[j] + 5, film_max)), method='bounded')
        keys = ('combined', 'delta', 'shape', 'level', 'reflectance')
        out.append(dict(name=name, film_nm=float(b.x), **dict(zip(keys, misfit(meas, refl, met, float(b.x), n_film, lam_range))),
                        bare=dict(zip(keys, misfit(meas, refl, met, 0.0, n_film, lam_range)))))
    return sorted(out, key=lambda x: x['combined'])


def report_identify(name, refl, ranking, say=print):
    lam = refl['channels']; t = int(np.argmin(np.abs(refl['theta_d'] - 30)))
    say(f'== {name}: luminance reflectance of the lobe at theta_d {refl["theta_d"][0]:.0f}-{refl["theta_d"][-1]:.0f} deg {np.mean(refl["level"]):.2f}; '
        'spectrum (luminance 1) at 462 / 550 / 646 / 798 nm: ' + ' / '.join(f'{refl["S"][t, int(np.argmin(np.abs(lam - x)))]:.2f}' for x in (462, 550, 646, 798)))
    say('  candidate                  film [nm]  combined  Delta [deg]  shape [log]  level [log]  reflectance | bare: combined  Delta')
    for r in ranking:
        say(f'  {r["name"]:26s} {r["film_nm"]:7.0f}  {r["combined"]:8.2f}  {r["delta"]:9.1f}  {r["shape"]:11.2f}  {r["level"]:11.2f}  {r["reflectance"]:11.2f} |'
            f' {r["bare"]["combined"]:13.2f}  {r["bare"]["delta"]:6.1f}')
