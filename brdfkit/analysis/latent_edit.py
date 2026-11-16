"""
Material editing in a latent space: move the grey angular shape A of a factor model (`brdfkit factor`) along the
inter-material PCA of MERL (`brdfkit merl-pca`), keeping the measured spectrum S and polarization N.

  1. A is put on the (strided) grid of the PCA model: phi_d is folded to [0, 180) deg by averaging the four equivalent
     angles phi, -phi, phi - 180, 180 - phi (isotropy, reciprocity), the theta axes are taken at the nearest bin.
     Grey -> RGB is (A, A, A).
  2. Its coefficients z are fitted (masked least squares, k components).
  3. An edit is a shift dz of the coefficients: `pcs` = [(p, t)] adds t * sigma_p along component p (1-based; sigma_p is the
     spread of the PCA's own materials along it), `toward` = [(name, t)] moves the fraction t of the way to a material of the PCA.
  4. Only the change is applied, in the log domain of the 'log' mapping:  A' = A * exp(mean_RGB(basis dz)).
     A itself is not replaced by its reconstruction, so what the basis cannot represent (the hpBRDF detail) is kept.
     The change is interpolated from the PCA grid to the grid of A (phi_d and theta_d linearly, theta_h linearly in
     sqrt(theta_h), the coordinate in which both grids are uniform).
  5. Optionally A' is rescaled so that its directional albedo at one angle of incidence equals that of A
     (`keep_albedo_deg`): only the shape changes, not the brightness.

The edited model is an ordinary factor model (.npz): `brdfkit compress` / `export` turn it into a .hpbpca for a renderer.
Large moves leave the space of plausible materials (the albedo can exceed 1); see docs/LATENT_EDIT.md.
"""
import warnings
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from .. import geometry as G
from .merl_pca import MaterialPCA

HALF_PI = np.pi / 2


def _sqrt_h(theta_h):
    return np.sqrt(np.clip(np.asarray(theta_h, float) / HALF_PI, 0, 1))


def _nearest(axis, x):
    return np.abs(np.asarray(axis)[None, :] - np.asarray(x)[:, None]).argmin(1)


def to_model_grid(model_f, pca: MaterialPCA):
    """grey A of the factor model on the strided grid of `pca` -> [P', T', H'] (NaN where A is unmeasured)"""
    A = np.asarray(model_f['A'], np.float64)
    phi = np.asarray(model_f['phi_d'], float)
    it = _nearest(model_f['theta_d'], pca.theta_d)
    ih = _nearest(_sqrt_h(model_f['theta_h']), _sqrt_h(pca.theta_h))
    out = np.full((len(pca.phi_d), len(it), len(ih)), np.nan)
    wrap = lambda a: (a + np.pi) % (2 * np.pi) - np.pi
    for j, p in enumerate(pca.phi_d):
        eq = [wrap(a) for a in (p, -p, p - np.pi, np.pi - p)]
        ip = sorted({int(np.abs(wrap(phi - a)).argmin()) for a in eq})
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)          # all-NaN bins (below the horizon) stay NaN
            out[j] = np.nanmean(np.stack([A[i][np.ix_(it, ih)] for i in ip]), 0)
    return out


def project(model_f, pca: MaterialPCA, k=16):
    """-> (coefficients [k], info) of the grey A on the PCA (masked least squares)"""
    g = to_model_grid(model_f, pca)
    sp, sd, sh = pca.stride
    v = np.full((g.shape[0] * sp, g.shape[1] * sd, g.shape[2] * sh, 3), np.nan)
    v[::sp, ::sd, ::sh] = g[..., None]
    return pca.project_masked(v, k)


def log_change(pca: MaterialPCA, dz):
    """change of log f on the strided PCA grid [P', T', H'] for a coefficient shift dz (grey: mean over RGB); bins outside the
    PCA mask take the last valid value along theta_h (0 where a whole row is outside)"""
    if pca.mapping != 'log':
        raise ValueError('latent editing needs a PCA with the log mapping')
    dz = np.asarray(dz, float)
    dx = (pca.basis[:, :len(dz)] @ dz).reshape(-1, 3).mean(1)
    D = np.zeros(pca.grid_shape); D[pca.valid] = dx
    for a in range(D.shape[0]):
        for b in range(D.shape[1]):
            v = pca.valid[a, b]
            if v.any() and not v.all():
                last = np.nonzero(v)[0][-1]; D[a, b, last + 1:] = D[a, b, last]
    return D


def to_factor_grid(D, pca: MaterialPCA, model_f):
    """interpolate a field on the strided PCA grid to the grid of the factor model [n_phi, n_thd, n_thh]"""
    gp = np.append(pca.phi_d, np.pi)                                  # phi_d = pi is phi_d = 0 (reciprocity)
    Dp = np.concatenate([D, D[:1]], 0)
    gd = np.asarray(pca.theta_d, float); gh = _sqrt_h(pca.theta_h)
    f = RegularGridInterpolator((gp, gd, gh), Dp, bounds_error=False, fill_value=None)
    phi = np.mod(np.abs(np.asarray(model_f['phi_d'], float)), np.pi)
    td = np.clip(np.asarray(model_f['theta_d'], float), gd[0], gd[-1]); th = np.clip(_sqrt_h(model_f['theta_h']), gh[0], gh[-1])
    P, T, H = np.meshgrid(phi, td, th, indexing='ij')
    return f(np.stack([P, T, H], -1))


def albedo(model_f, theta_i_deg, n=400_000, seed=0):
    """directional albedo (luminance, since S is normalized to Y = 1) of the grey A at one angle of incidence: Monte Carlo over
    the hemisphere with nearest-bin lookups"""
    A = np.asarray(model_f['A'], np.float64)
    rng = np.random.default_rng(seed); u1, u2 = rng.random(n), rng.random(n)
    ct = u1; st = np.sqrt(1 - ct ** 2); ph = 2 * np.pi * u2
    o = np.stack([st * np.cos(ph), st * np.sin(ph), ct], -1); t = np.radians(theta_i_deg)
    i = np.ascontiguousarray(np.broadcast_to([np.sin(t), 0.0, np.cos(t)], o.shape))
    pd, th, td = G.rusinkiewicz(i, o)
    ip = _nearest(model_f['phi_d'], pd); it = _nearest(model_f['theta_d'], td); ih = _nearest(_sqrt_h(model_f['theta_h']), _sqrt_h(th))
    f = A[ip, it, ih]
    return float(np.mean(np.where(np.isfinite(f), f, 0.0) * o[:, 2]) * 2 * np.pi)


def edit(model_f, pca: MaterialPCA, pcs=(), toward=(), k=16, keep_albedo_deg=None):
    """-> (edited factor model, info). pcs: [(component 1-based, t in units of the PCA's spread)], toward: [(material, fraction)]"""
    k = min(int(k), pca.k)
    z, pinfo = project(model_f, pca, k)
    sig = pca.coef[:, :k].std(0)
    dz = np.zeros(k)
    for p, t in pcs:
        if not 1 <= int(p) <= k:
            raise ValueError(f'component {p} outside 1..{k}')
        dz[int(p) - 1] += float(t) * sig[int(p) - 1]
    for name, t in toward:
        if name not in pca.names:
            raise ValueError(f'{name!r} is not a material of the PCA')
        dz += float(t) * (pca.coef[pca.names.index(name), :k] - z)
    delta = to_factor_grid(log_change(pca, dz), pca, model_f)
    A = np.asarray(model_f['A'], np.float64)
    A2 = A * np.exp(delta)
    out = dict(model_f); out['A'] = A2.astype(np.float32)
    scale = 1.0
    if keep_albedo_deg is not None:
        scale = albedo(model_f, keep_albedo_deg) / max(albedo(out, keep_albedo_deg), 1e-300)
        out['A'] = (A2 * scale).astype(np.float32)
    z2 = z + dz
    near = [pca.names[i] for i in np.argsort(np.linalg.norm(pca.coef[:, :k] - z2, axis=1))[:5]]
    info = dict(k=k, z=(z / sig).tolist(), dz=(dz / sig).tolist(), projection_rmse=pinfo['rmse'], scale=scale, nearest=near,
                albedo_before=[albedo(model_f, a) for a in (0, 45, 75)], albedo_after=[albedo(out, a) for a in (0, 45, 75)],
                pcs=[[int(p), float(t)] for p, t in pcs], toward=[[n, float(t)] for n, t in toward], keep_albedo_deg=keep_albedo_deg)
    out['meta'] = dict(model_f.get('meta', {}), edit=info)
    return out, info


def report(info, say=print):
    say(f'coefficients of A on the PCA (in units of the spread of its materials, k = {info["k"]}, rmse {info["projection_rmse"]:.3f}):')
    say('  ' + ' '.join(f'{v:+.1f}' for v in info['z']))
    say('shift: ' + ', '.join([f'PC{p} {t:+g}' for p, t in info['pcs']] + [f'{t:g} of the way to {n}' for n, t in info['toward']]))
    say(f'albedo (luminance) at 0 / 45 / 75 deg: before ' + ' '.join(f'{a:.3f}' for a in info['albedo_before']) +
        '   after ' + ' '.join(f'{a:.3f}' for a in info['albedo_after']) +
        (f'   (rescaled by {info["scale"]:.2f} to keep the albedo at {info["keep_albedo_deg"]:g} deg)' if info['keep_albedo_deg'] is not None else ''))
    if max(info['albedo_after']) > 1:
        say('  warning: albedo above 1, the edited material creates energy; use a smaller shift or --keep-albedo')
    say('nearest materials of the PCA after the shift: ' + ', '.join(info['nearest']))
