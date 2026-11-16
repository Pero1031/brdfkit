"""
Factorized model of a polarized spectral table (hpBRDF):

    m00(phi_d, theta_d, theta_h, lambda) = A(phi_d, theta_d, theta_h) * S(theta_d, lambda)
    m_ij / m00                           = N_ij(theta_d, lambda)                     (15 elements)

A is a grey angular shape in luminance units (S is normalized to luminance Y = 1 for every theta_d), so it can be compared
with or replaced by a MERL-like lobe; S is the spectrum and N the polarization, both functions of theta_d (the angle of
incidence on the microfacet) and wavelength only. The table's Mueller frames have their x reference vectors in the plane of
reflection, which for a mirror microfacet is its plane of incidence, so a single-bounce specular reflection has exactly this
form (Fresnel in (theta_d, lambda)). See docs/FACTOR.md for how well the hpBRDF tables follow it.

Fit (two streaming passes over the table):
  pass 1  per theta_d: Gram matrix G = sum over bins of m00 m00^T (n_ch x n_ch) -> S = its top eigenvector (the energy-optimal
          rank-1 spectrum); N = m00-weighted mean of clip(m_ij / m00, -1, 1) over phi_d and theta_h.
          Outside the measured theta_d range (hpBRDF: about 20-80 deg, the rest is inpainted) S and N are copied from the
          nearest measured theta_d. S may be smoothed along theta_d (Gaussian, `smooth_deg`), then normalized to Y = 1.
  pass 2  A = m00 . S / (S . S) per bin, and the reconstruction error of m00 (energy weighted) and of N, by theta_h region.
Several S variants (smoothing widths) are fitted in the same two passes.

Model file (.npz): kind = 'brdfkit-factor', A [n_phi, n_thd, n_thh] float16, S [n_thd, n_ch], N [n_thd, n_ch, 15],
phi_d, theta_d, theta_h, channels, measured_theta_d, smooth_deg, errors (json), meta (json).
`view(model)` turns a model back into a BRDFTable (slices built on access), e.g. for `brdfkit compress`.
"""
import json
import os
import time
import numpy as np
from scipy.ndimage import gaussian_filter1d
from ..table import BRDFTable, WAVELENGTH
from ..color import cmf

KIND = 'brdfkit-factor'
REGIONS = ((0.0, 5.0), (5.0, 20.0), (20.0, 90.01))      # theta_h ranges [deg] for the error report
RIDX = np.arange(16)[1:]                                 # Mueller elements other than m00, row-major


def luminance_weights(channels):
    """band weights giving the CIE Y of a spectrum sampled at `channels` (sum to 1, so a flat spectrum has Y = its value)"""
    y = cmf(channels)[:, 1]
    return y / y.sum()


def _region_index(theta_h):
    th = np.degrees(theta_h); r = np.full(len(th), len(REGIONS) - 1)
    for k, (a, b) in enumerate(REGIONS):
        r[(th >= a) & (th < b)] = k
    return r


def _extend_edges(X, measured):
    """copy the first / last measured theta_d row of X [n_thd, ...] into the unmeasured rows outside the measured range"""
    idx = np.nonzero(measured)[0]
    X = X.copy(); X[:idx[0]] = X[idx[0]]; X[idx[-1] + 1:] = X[idx[-1]]
    return X


def fit(table: BRDFTable, measured_theta_d=(20.0, 80.0), smooth_deg=(0.0,), log=print):
    """-> list of model dicts, one per smoothing width in `smooth_deg` (degrees of theta_d; 0 = none)"""
    if not table.polarized or table.channel_kind != WAVELENGTH:
        raise ValueError('factor needs a polarized spectral table')
    nP, nT, nH, nC = table.n_phi, table.n_thd, table.n_thh, table.n_ch
    tdd = np.degrees(table.theta_d)
    meas = (tdd >= measured_theta_d[0] - 1e-6) & (tdd <= measured_theta_d[1] + 1e-6)
    if not meas.any():
        raise ValueError('no theta_d inside the measured range')
    reg = _region_index(table.theta_h)
    G = np.zeros((nT, nC, nC)); Wn = np.zeros((nT, nC)); Sn = np.zeros((nT, nC, 15))
    t0 = time.time()
    for i, sl in table.slices():
        m = sl[..., 0, 0].astype(np.float64)                                   # [T, H, C]
        ok = np.isfinite(sl.reshape(nT, nH, -1)).all(-1) & (m > 0).all(-1)
        mz = np.where(ok[..., None], m, 0.0)
        G += np.einsum('thc,thd->tcd', mz, mz)
        with np.errstate(invalid='ignore', divide='ignore'):
            N = np.clip(sl.reshape(nT, nH, nC, 16)[..., RIDX] / m[..., None], -1, 1)
        N = np.where(ok[..., None, None], N, 0.0)
        Wn += mz.sum(1); Sn += (mz[..., None] * N).sum(1)
        if i % 60 == 0:
            log(f'  pass 1  {i}/{nP}  {time.time() - t0:.0f}s')
    S_raw = np.zeros((nT, nC))
    for j in range(nT):
        if G[j].trace() > 0:
            w, V = np.linalg.eigh(G[j]); v = V[:, -1]; S_raw[j] = v * np.sign(v.sum())
    S_raw = _extend_edges(S_raw, meas)
    N_fit = _extend_edges(Sn / np.maximum(Wn, 1e-300)[..., None], meas)
    yw = luminance_weights(table.channels)
    variants = []
    for sd in smooth_deg:
        S = S_raw.copy()
        if sd > 0:
            step = np.median(np.diff(tdd)); idx = np.nonzero(meas)[0]
            S[idx] = gaussian_filter1d(S[idx], sd / step, axis=0, mode='nearest')
            S = _extend_edges(S, meas)
        S = S / (S @ yw)[:, None]
        variants.append(dict(smooth=float(sd), S=S, A=np.full((nP, nT, nH), np.nan, np.float32)))
    log(f'pass 1 done {time.time() - t0:.0f}s')

    # pass 2: A and the errors
    nR = len(REGIONS)
    e = [dict(num=np.zeros((2, nR)), den=np.zeros((2, nR))) for _ in variants]      # [all / measured theta_d, region]
    en = dict(num=np.zeros((2, nR)), den=np.zeros((2, nR)))                          # N error, energy weighted
    t0 = time.time()
    for i, sl in table.slices():
        m = sl[..., 0, 0].astype(np.float64)
        ok = np.isfinite(sl.reshape(nT, nH, -1)).all(-1) & (m > 0).all(-1)
        mz = np.where(ok[..., None], m, 0.0)
        e2 = mz ** 2
        for v, ev in zip(variants, e):
            S = v['S']
            A = np.einsum('thc,tc->th', mz, S) / (S * S).sum(1)[:, None]
            v['A'][i] = np.where(ok, A, np.nan)
            r2 = ((mz - A[..., None] * S[:, None, :]) ** 2)
            for k in range(nR):
                sel = ok & (reg[None, :] == k)
                for a, rows in enumerate((np.ones(nT, bool), meas)):
                    s2 = sel & rows[:, None]
                    ev['num'][a, k] += r2[s2].sum(); ev['den'][a, k] += e2[s2].sum()
        with np.errstate(invalid='ignore', divide='ignore'):
            N = np.clip(sl.reshape(nT, nH, nC, 16)[..., RIDX] / m[..., None], -1, 1)
        dN = np.where(ok[..., None, None], (N - N_fit[:, None]) ** 2, 0.0).sum(-1) * mz   # [T, H, C]
        for k in range(nR):
            sel = ok & (reg[None, :] == k)
            for a, rows in enumerate((np.ones(nT, bool), meas)):
                s2 = sel & rows[:, None]
                en['num'][a, k] += dN[s2].sum(); en['den'][a, k] += 15 * mz[s2].sum()
        if i % 60 == 0:
            log(f'  pass 2  {i}/{nP}  {time.time() - t0:.0f}s')
    log(f'pass 2 done {time.time() - t0:.0f}s')

    def rel(d, a):
        return dict(**{f'theta_h {lo:g}-{min(hi, 90):g}': float(np.sqrt(d['num'][a, k] / max(d['den'][a, k], 1e-300)))
                       for k, (lo, hi) in enumerate(REGIONS)}, all=float(np.sqrt(d['num'][a].sum() / max(d['den'][a].sum(), 1e-300))))
    models = []
    for v, ev in zip(variants, e):
        errors = dict(m00_rel_rmse=dict(all_theta_d=rel(ev, 0), measured_theta_d=rel(ev, 1)),
                      N_rmse=dict(all_theta_d=rel(en, 0), measured_theta_d=rel(en, 1)))
        models.append(dict(kind=KIND, A=v['A'], S=v['S'], N=N_fit, phi_d=table.phi_d, theta_d=table.theta_d, theta_h=table.theta_h,
                           channels=table.channels, measured_theta_d=np.array(measured_theta_d, float), smooth_deg=v['smooth'],
                           errors=errors, meta=dict(source=str(table.meta.get('source', '')),
                                                   floor_offset=table.meta.get('floor_offset'), floor_max_gain=table.meta.get('floor_max_gain'))))
    return models


def save(path, model):
    z = dict(model)
    z['A'] = model['A'].astype(np.float16)
    z['errors'] = json.dumps(model['errors']); z['meta'] = json.dumps(model['meta'], default=str)
    np.savez(path, **z)
    return path


def load(path):
    with np.load(path, allow_pickle=False) as z:
        m = {k: z[k] for k in z.files}
    if str(m.get('kind')) != KIND:
        raise ValueError(f'{path}: not a brdfkit factor model')
    m['errors'] = json.loads(str(m['errors'])); m['meta'] = json.loads(str(m['meta']))
    m['smooth_deg'] = float(m['smooth_deg'])
    return m


def is_model(path):
    if not str(path).endswith('.npz') or not os.path.exists(path):
        return False
    with np.load(path, allow_pickle=False) as z:
        return 'kind' in z.files and str(z['kind']) == KIND


def view(model) -> BRDFTable:
    """BRDFTable of the model (polarized, spectral); each phi_d slice is rebuilt on access"""
    if isinstance(model, str):
        model = load(model)
    A, S, N = model['A'], np.asarray(model['S'], np.float32), np.asarray(model['N'], np.float32)
    nT, nC = S.shape

    def slice_fn(i):
        m00 = A[i].astype(np.float32)[..., None] * S[:, None, :]                  # [T, H, C]
        M = np.empty(m00.shape + (16,), np.float32)
        M[..., 0] = m00; M[..., 1:] = m00[..., None] * N[:, None]
        return M.reshape(m00.shape + (4, 4))
    meta = dict(source=f'factor model of {model["meta"].get("source", "")}', smooth_deg=float(model['smooth_deg']))
    return BRDFTable(model['phi_d'], model['theta_d'], model['theta_h'], model['channels'], WAVELENGTH, True, meta, slice_fn=slice_fn)


def report(model, say=print):
    say(f'factor model of {model["meta"].get("source", "")}  (S smoothed over theta_d: {model["smooth_deg"]:g} deg; '
        f'measured theta_d {model["measured_theta_d"][0]:g}-{model["measured_theta_d"][1]:g} deg)')
    for key, lab in (('m00_rel_rmse', 'm00 relative RMS error (energy weighted)'), ('N_rmse', 'N = m_ij/m00 RMS error (m00 weighted)')):
        for scope in ('measured_theta_d', 'all_theta_d'):
            d = model['errors'][key][scope]
            say(f'  {lab:42s} {scope:17s} ' + '  '.join(f'{k} {v:.4f}' for k, v in d.items()))
