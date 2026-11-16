"""
Inter-material PCA of a set of BRDF tables on one common grid (Matusik, Pfister, Brand, McMillan 2003, section 4).

Each material becomes one feature vector x over the measured bins of the (strided) grid times the 3 colour channels:
    mapping 'log'     : x = log(max(f, eps))                                   (Matusik 2003)
    mapping 'nielsen' : x = log((f cos_i cos_o + eps) / (median_m f cos_i cos_o + eps))   (Nielsen, Jensen, Ramamoorthi 2015)
The PCA is solved through the n x n Gram matrix of the centred vectors (n = number of materials), which gives every eigenvalue,
the coefficients of every material, and the leave-one-out residual of each material against the basis fitted on the others,
without ever forming a D x D covariance (D ~ 5e5 .. 4e6).

Result (`MaterialPCA`, saved as .npz): names, mean [D], basis [D, k], coef [n, n] (all components), eigvals, evr, valid mask
[P', T', H'], the strided axes, mapping, eps, ref [D] (nielsen median) and cosw [P', T', H'] (cos_i cos_o).
`project(values)` maps a table on the same (unstrided) grid to coefficients; `reconstruct(coef)` returns linear BRDF values on the
strided grid. `run_report` fits a directory of MERL .binary files and writes the npz, a text report and figures.
"""
import glob
import json
import os
import time
import numpy as np
from ..table import BRDFTable

EPS = {'log': 1e-5, 'nielsen': 1e-3}
CHUNK = 1 << 16


# ----------------------------------------------------------------------------- grid
def _stride3(stride):
    s = (stride,) * 3 if isinstance(stride, int) else tuple(int(v) for v in stride)
    if len(s) != 3 or min(s) < 1:
        raise ValueError('stride must be an int or 3 ints (phi_d, theta_d, theta_h)')
    return s


def cos_weights(phi_d, theta_d, theta_h):
    """cos(theta_i) cos(theta_o) on the grid [P, T, H]; 0 where either direction is below the horizon
    (MERL half/diff convention: d rotated by theta_h about y gives omega_i, omega_o = 2 (omega_i . h) h - omega_i)"""
    P, T, H = np.meshgrid(phi_d, theta_d, theta_h, indexing='ij')
    dx, dz = np.sin(T) * np.cos(P), np.cos(T)
    ci = -dx * np.sin(H) + dz * np.cos(H)                                     # omega_i . z
    co = 2 * np.cos(T) * np.cos(H) - ci                                       # omega_o . z   (omega_i . h = cos theta_d)
    w = ci * co
    w[(ci <= 0) | (co <= 0)] = 0.0
    return w


# ----------------------------------------------------------------------------- the model
class MaterialPCA:
    KEYS = ('names', 'mean', 'basis', 'coef', 'eigvals', 'evr', 'valid', 'phi_d', 'theta_d', 'theta_h', 'stride', 'mapping',
            'eps', 'ref', 'cosw', 'recon_rmse', 'loo_rmse', 'meta')

    def __init__(self, **kw):
        for k in self.KEYS:
            setattr(self, k, kw.get(k))
        self.names = list(self.names); self.mapping = str(self.mapping); self.eps = float(self.eps)
        self.stride = _stride3(tuple(int(v) for v in np.atleast_1d(self.stride)))
        if isinstance(self.meta, (str, np.ndarray)):
            self.meta = json.loads(str(self.meta))

    # ---- sizes
    @property
    def n(self): return len(self.names)
    @property
    def k(self): return self.basis.shape[1]
    @property
    def dim(self): return self.basis.shape[0]
    @property
    def grid_shape(self): return self.valid.shape

    # ---- features
    def _subsample(self, values):
        sp, sd, sh = self.stride
        v = np.asarray(values)[::sp, ::sd, ::sh]
        if v.shape[:3] != self.grid_shape or v.shape[-1] != 3:
            raise ValueError(f'table grid {v.shape} does not match the model grid {self.grid_shape} + (3,) after stride {self.stride}')
        return v

    def features(self, values):
        """table values [P, T, H, 3] on the unstrided grid -> feature vector [D] (float64); NaN where the table is unmeasured"""
        v = self._subsample(values)[self.valid].astype(np.float64)               # [n_valid, 3]
        if self.mapping == 'log':
            x = np.log(np.maximum(v, self.eps))
        else:
            w = self.cosw[self.valid][:, None]
            x = np.log((v * w + self.eps) / (self.ref.reshape(-1, 3) + self.eps))
        x[~np.isfinite(v)] = np.nan
        return x.ravel()

    def unmap(self, x):
        """feature vector [D] -> linear BRDF values [n_valid, 3]"""
        x = x.reshape(-1, 3)
        if self.mapping == 'log':
            return np.exp(x)
        w = self.cosw[self.valid][:, None]
        return np.maximum((np.exp(x) * (self.ref.reshape(-1, 3) + self.eps) - self.eps), 0.0) / np.where(w > 0, w, np.inf)

    def project(self, values, k=None):
        """coefficients [k] of one table (values [P, T, H, 3]; unmeasured bins inside the common mask are not allowed)"""
        x = self.features(values)
        if not np.isfinite(x).all():
            raise ValueError('table has unmeasured bins inside the model mask')
        return (x - self.mean) @ self.basis[:, :k]

    def project_masked(self, values, k=None):
        """least-squares coefficients [k] of a table that may be unmeasured (NaN) in part of the model mask
        (the full-mask `project` is the special case where the basis rows are orthonormal).
        -> (coef, info) with info = n_used, n_total, cond (condition number of the restricted basis), rmse (feature-space RMSE over the used entries)"""
        k = self.k if k is None else min(int(k), self.k)
        x = self.features(values) - self.mean
        use = np.isfinite(x)
        if use.sum() < 8 * k:
            raise ValueError(f'only {use.sum()} measured features inside the model mask, too few for {k} coefficients')
        B = self.basis[use, :k]
        coef, _, _, sv = np.linalg.lstsq(B, x[use], rcond=None)
        rmse = float(np.sqrt(np.mean((x[use] - B @ coef) ** 2)))
        return coef, dict(n_used=int(use.sum()), n_total=int(use.size), cond=float(sv[0] / max(sv[-1], 1e-300)), rmse=rmse)

    def reconstruct(self, coef, k=None):
        """coefficients [k'] -> linear BRDF values on the strided grid [P', T', H', 3] (NaN outside the common mask)"""
        c = np.asarray(coef, np.float64)[:k]
        x = self.mean + self.basis[:, :len(c)] @ c
        out = np.full(self.grid_shape + (3,), np.nan, np.float32)
        out[self.valid] = self.unmap(x)
        return out

    def table(self, coef, k=None, meta=None) -> BRDFTable:
        from ..table import RGB
        return BRDFTable(self.phi_d, self.theta_d, self.theta_h, np.arange(3.0), RGB, False, meta or dict(source='merl_pca'),
                         values=self.reconstruct(coef, k))

    # ---- io
    def save(self, path):
        d = {k: getattr(self, k) for k in self.KEYS}
        d['names'] = np.array(self.names); d['stride'] = np.array(self.stride); d['meta'] = json.dumps(self.meta or {})
        d['ref'] = np.zeros(0) if self.ref is None else self.ref
        for k in ('mean', 'basis', 'ref', 'cosw', 'recon_rmse', 'loo_rmse'):
            d[k] = np.asarray(d[k], np.float32)
        np.savez(path, **d)
        return path

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as zf:
            d = {k: zf[k] for k in cls.KEYS if k in zf}
        d['names'] = [str(s) for s in d['names']]; d['mapping'] = str(d['mapping'])
        if d.get('ref') is not None and d['ref'].size == 0:
            d['ref'] = None
        for k in ('mean', 'basis', 'coef', 'eigvals', 'evr'):
            d[k] = d[k].astype(np.float64)
        return cls(**d)


# ----------------------------------------------------------------------------- fitting
def list_materials(paths):
    """directory or list of files -> list of (name, path), sorted by name"""
    if isinstance(paths, str):
        paths = [paths]
    files = []
    for p in paths:
        files += sorted(glob.glob(os.path.join(p, '*.binary'))) if os.path.isdir(p) else [p]
    return sorted((os.path.splitext(os.path.basename(f))[0], f) for f in files)


def _gram(X):
    """X [n, D] float32 -> X X^T in float64, accumulated in column chunks"""
    n = len(X); G = np.zeros((n, n))
    for j in range(0, X.shape[1], CHUNK):
        B = X[:, j:j + CHUNK].astype(np.float64)
        G += B @ B.T
    return G


def _loo(G, n_comp):
    """leave-one-out squared residual of each material against the PCA of the others, for every k (from the centred Gram only)
    -> resid2 [n, n_comp + 1] (column k = residual after k components; column 0 = ||y||^2)"""
    n = len(G); out = np.zeros((n, n_comp + 1))
    for i in range(n):
        J = np.array([j for j in range(n) if j != i]); m = len(J)
        GJ = G[np.ix_(J, J)]
        rowm = GJ.mean(1); allm = GJ.mean()
        GJc = GJ - rowm[:, None] - rowm[None, :] + allm
        w, V = np.linalg.eigh(GJc); order = np.argsort(w)[::-1]; w = np.maximum(w[order], 0); V = V[:, order]
        b = G[i, J] - G[i, J].mean() - rowm + allm
        y2 = G[i, i] - 2 * G[i, J].mean() + allm
        c = (b @ V[:, :n_comp]) / np.sqrt(np.maximum(w[:n_comp], 1e-300))
        c[w[:n_comp] <= 1e-12 * max(w[0], 1e-300)] = 0.0
        out[i, 0] = y2; out[i, 1:] = y2 - np.cumsum(c ** 2)
    return np.maximum(out, 0)


def fit(materials, stride=2, mapping='log', eps=None, k=16, loo=True, log=print):
    """materials: list of (name, path) or (name, BRDFTable). Returns (MaterialPCA, Xc) where Xc [n, D] float32 are the centred
    feature vectors (kept for plots; discard if memory matters)."""
    from ..io import load as load_table
    if mapping not in EPS:
        raise ValueError("mapping must be 'log' or 'nielsen'")
    eps = EPS[mapping] if eps is None else float(eps)
    stride = _stride3(stride); sp, sd, sh = stride
    names = [m[0] for m in materials]; n = len(names)
    t0 = time.time(); X = None; valid = None; axes = None
    for i, (name, src) in enumerate(materials):
        t = src if isinstance(src, BRDFTable) else load_table(src)
        if t.polarized or t.n_ch != 3:
            raise ValueError(f'{name}: inter-material PCA expects unpolarized RGB tables')
        v = t.values[::sp, ::sd, ::sh]
        if X is None:
            axes = (t.phi_d[::sp], t.theta_d[::sd], t.theta_h[::sh])
            X = np.empty((n,) + v.shape, np.float32); valid = np.ones(v.shape[:3], bool)
        elif v.shape != X.shape[1:]:
            raise ValueError(f'{name}: grid {v.shape} differs from {X.shape[1:]}')
        X[i] = v; valid &= np.isfinite(v).all(-1)
        if i % 20 == 0 or i == n - 1:
            log(f'loaded {i + 1}/{n}  {time.time() - t0:.0f}s')
    cosw = cos_weights(*axes)
    if mapping == 'nielsen':
        valid &= cosw > 0
    log(f'grid {valid.shape} (stride {stride}), common measured bins {valid.mean():.3f} -> D = {valid.sum() * 3}')
    Xv = np.ascontiguousarray(X[:, valid].reshape(n, -1)); del X                      # [n, D] linear
    ref = None
    if mapping == 'log':
        np.log(np.maximum(Xv, eps, out=Xv), out=Xv)
    else:
        w = np.repeat(cosw[valid], 3).astype(np.float32)
        Xv *= w
        ref = np.median(Xv, axis=0)
        Xv += eps; Xv /= (ref + eps); np.log(Xv, out=Xv)
    mean = Xv.mean(0, dtype=np.float64)
    Xv -= mean.astype(np.float32)
    D = Xv.shape[1]
    t0 = time.time(); G = _gram(Xv)
    w, V = np.linalg.eigh(G); order = np.argsort(w)[::-1]; w = np.maximum(w[order], 0); V = V[:, order]
    evr = w / max(w.sum(), 1e-300)
    sw = np.sqrt(w)
    coef = V * sw                                                                     # [n, n] coefficients of every component
    k = min(k, n - 2)                                                                 # n - 2: leaves a LOO column for k
    basis = np.zeros((D, k))
    for j in range(0, D, CHUNK):
        basis[j:j + CHUNK] = Xv[:, j:j + CHUNK].astype(np.float64).T @ (V[:, :k] / np.maximum(sw[:k], 1e-300))
    recon2 = np.diag(G)[:, None] - np.concatenate([np.zeros((n, 1)), np.cumsum(coef ** 2, 1)], 1)   # [n, n+1]
    recon_rmse = np.sqrt(np.maximum(recon2, 0) / D)
    loo_rmse = np.sqrt(_loo(G, n - 2) / D) if loo and n >= 3 else np.zeros((n, 0))
    log(f'PCA (Gram {n}x{n}) {time.time() - t0:.1f}s;  cum. explained variance: ' +
        ', '.join(f'k={kk}: {evr[:kk].sum():.4f}' for kk in (1, 2, 3, 5, 8, 16, 32) if kk < n))
    meta = dict(n_materials=n, dim=D, stride=list(stride), mapping=mapping, eps=eps, sources=[m[1] if isinstance(m[1], str) else 'table' for m in materials],
                features='x = log(max(f, eps))' if mapping == 'log' else 'x = log((f cos_i cos_o + eps) / (median_m f cos_i cos_o + eps))',
                layout='mean/basis/ref over bins valid[P,T,H] (C order) x 3 channels; coef[n, n] all components; basis[D, k]')
    model = MaterialPCA(names=names, mean=mean, basis=basis, coef=coef, eigvals=w, evr=evr, valid=valid, phi_d=axes[0], theta_d=axes[1],
                        theta_h=axes[2], stride=stride, mapping=mapping, eps=eps, ref=ref, cosw=cosw, recon_rmse=recon_rmse,
                        loo_rmse=loo_rmse, meta=meta)
    return model, Xv


# ----------------------------------------------------------------------------- report
def linear_error(model, Xc, ks):
    """relative RMSE of cos-weighted linear BRDF (f cos_i cos_o) per material for each k in ks -> [n, len(ks)]
    (reconstruction of the fitted material; measures what the mapping and truncation cost in radiance terms)"""
    w = np.repeat(model.cosw[model.valid], 3)
    out = np.zeros((model.n, len(ks)))
    for i in range(model.n):
        x = Xc[i].astype(np.float64)
        f = (model.unmap(model.mean + x).ravel() * w)
        nrm = np.sqrt(np.mean(f ** 2))
        for j, k in enumerate(ks):
            r = model.unmap(model.mean + model.basis[:, :k] @ model.coef[i, :k]).ravel() * w
            out[i, j] = np.sqrt(np.mean((r - f) ** 2)) / nrm
    return out


def run_report(paths, out_dir, stride=2, k=16, mapping='log', eps=None, plot=True, log=print):
    """fit every .binary under `paths`, write merl_pca.npz, merl_pca_report.txt and figures to out_dir"""
    materials = list_materials(paths)
    if len(materials) < 3:
        raise ValueError(f'need at least 3 materials, found {len(materials)}')
    os.makedirs(out_dir, exist_ok=True)
    model, Xc = fit(materials, stride, mapping, eps, k, log=log)
    npz = model.save(os.path.join(out_dir, 'merl_pca.npz'))
    n, K = model.n, model.k
    report = []
    def say(*a):
        s = ' '.join(str(x) for x in a); log(s); report.append(s)

    say(f'=== inter-material PCA: {n} materials, mapping {mapping} (eps {model.eps:g}), stride {model.stride}, D = {model.dim} ===')
    say(f'sources: {paths}')
    ks = [kk for kk in (1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64) if kk <= n - 2]
    say('\nk   cum.var   recon RMSE (feature space) median / p90 / max   LOO RMSE median / p90 / max')
    for kk in ks:
        r = model.recon_rmse[:, kk]; l = model.loo_rmse[:, kk]
        say(f'{kk:3d}  {model.evr[:kk].sum():.4f}   {np.median(r):.3f} / {np.percentile(r, 90):.3f} / {r.max():.3f}'
            f'                      {np.median(l):.3f} / {np.percentile(l, 90):.3f} / {l.max():.3f}')
    say(f'(feature-space RMSE is in natural-log units: 0.1 = 10% per-bin error; LOO = material held out of the fit)')
    lk = [kk for kk in (1, 2, 4, 8, K) if kk <= K]
    le = linear_error(model, Xc, lk)
    say('\nrelative RMSE of f cos_i cos_o (linear, cos-weighted), reconstruction of the fitted material:')
    say('k    ' + '  '.join(f'{kk:>6d}' for kk in lk))
    say('med  ' + '  '.join(f'{np.median(le[:, j]):6.3f}' for j in range(len(lk))))
    say('p90  ' + '  '.join(f'{np.percentile(le[:, j], 90):6.3f}' for j in range(len(lk))))
    say('max  ' + '  '.join(f'{le[:, j].max():6.3f}' for j in range(len(lk))))

    say(f'\nper material: first 4 coefficients, recon RMSE at k={K}, LOO RMSE at k={K}, nearest neighbour in the {K}-d space')
    C = model.coef[:, :K]; dist = np.sqrt(((C[:, None] - C[None]) ** 2).sum(-1)); np.fill_diagonal(dist, np.inf)
    nn = dist.argmin(1)
    say(f'{"material":28s} {"c1":>8s} {"c2":>8s} {"c3":>8s} {"c4":>8s}  {"recon":>6s}  {"loo":>6s}  nearest')
    for i in np.argsort(model.loo_rmse[:, K]):
        say(f'{model.names[i]:28s} ' + ' '.join(f'{model.coef[i, j]:8.1f}' for j in range(4)) +
            f'  {model.recon_rmse[i, K]:6.3f}  {model.loo_rmse[i, K]:6.3f}  {model.names[nn[i]]} ({dist[i, nn[i]]:.1f})')
    say(f'\neigenvalues (variance per component, first 10): ' + ' '.join(f'{v:.4g}' for v in model.eigvals[:10]))
    open(os.path.join(out_dir, 'merl_pca_report.txt'), 'w', encoding='utf-8').write('\n'.join(report) + '\n')
    if plot:
        try:
            _plots(model, Xc, out_dir, le, lk)
        except ImportError as e:
            log('plots skipped:', e)
    log('wrote', npz)
    return model


def _slice_image(model, x, ipd=None):
    """feature vector [D] -> channel-mean image [T', H'] at one phi_d slice (NaN outside the mask)"""
    ipd = model.grid_shape[0] // 2 if ipd is None else ipd
    full = np.full(model.grid_shape + (3,), np.nan); full[model.valid] = x.reshape(-1, 3)
    return full[ipd].mean(-1)


def _plots(model, Xc, out_dir, le, lk):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n, K = model.n, model.k
    kk = np.arange(1, n - 1)

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.4))
    axs[0].semilogy(np.arange(1, n), np.maximum(model.evr[:n - 1], 1e-12), 'o-', ms=3, label='per component')   # n-1: centring removes one
    axs[0].set(title='explained variance ratio', xlabel='component', ylabel='fraction'); axs[0].grid(alpha=0.3)
    ax2 = axs[0].twinx(); ax2.plot(np.arange(1, n), np.cumsum(model.evr)[:n - 1], 'r-', label='cumulative'); ax2.set_ylim(0, 1.001); ax2.set_ylabel('cumulative', color='r')
    for q, ls, lab in ((50, '-', 'median'), (90, '--', 'p90'), (100, ':', 'max')):
        axs[1].plot(kk, np.percentile(model.recon_rmse[:, 1:n - 1], q, axis=0), 'b' + ls, label=f'reconstruction {lab}')
        if model.loo_rmse.shape[1] > 1:
            axs[1].plot(kk, np.percentile(model.loo_rmse[:, 1:n - 1], q, axis=0), 'r' + ls, label=f'leave-one-out {lab}')
    axs[1].set(title=f'RMSE in feature space ({model.mapping} mapping)', xlabel='k', ylabel='RMSE (ln units)', xscale='log'); axs[1].grid(alpha=0.3); axs[1].legend(fontsize=8)
    axs[2].plot(lk, np.median(le, 0), 'o-', label='median'); axs[2].plot(lk, np.percentile(le, 90, 0), 's--', label='p90'); axs[2].plot(lk, le.max(0), '^:', label='max')
    axs[2].set(title='relative RMSE of f cos_i cos_o (reconstruction)', xlabel='k', ylabel='rel. RMSE', xscale='log', yscale='log'); axs[2].grid(alpha=0.3); axs[2].legend(fontsize=8)
    fig.suptitle(f'inter-material PCA, {n} materials, D = {model.dim}'); fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'merl_pca_variance.png'), dpi=130); plt.close(fig)

    peak = np.array([np.log10(np.nanmax(model.unmap(model.mean + Xc[i])) / np.nanmedian(model.unmap(model.mean + Xc[i]))) for i in range(n)])
    fig, axs = plt.subplots(1, 2, figsize=(17, 8))
    for a, (i, j) in zip(axs, ((0, 1), (2, 3))):
        if j >= K: continue
        sc = a.scatter(model.coef[:, i], model.coef[:, j], c=peak, cmap='viridis', s=28, edgecolor='k', lw=0.3)
        for m in range(n):
            a.annotate(model.names[m], (model.coef[m, i], model.coef[m, j]), fontsize=6, xytext=(2, 2), textcoords='offset points')
        a.set(xlabel=f'PC{i + 1} ({model.evr[i]:.1%})', ylabel=f'PC{j + 1} ({model.evr[j]:.1%})'); a.grid(alpha=0.3); a.axhline(0, color='gray', lw=0.5); a.axvline(0, color='gray', lw=0.5)
    fig.colorbar(sc, ax=axs, shrink=0.6, label='log10(peak / median BRDF)  (specularity)')
    fig.suptitle('materials in the principal-component space'); fig.savefig(os.path.join(out_dir, 'merl_pca_scatter.png'), dpi=130); plt.close(fig)

    nshow = min(6, K); ipd = model.grid_shape[0] // 2
    ext = [np.degrees(model.theta_h[0]), np.degrees(model.theta_h[-1]), np.degrees(model.theta_d[0]), np.degrees(model.theta_d[-1])]
    fig, axs = plt.subplots(1, nshow + 1, figsize=(3.2 * (nshow + 1), 3.6))
    im = axs[0].imshow(_slice_image(model, model.mean, ipd), origin='lower', extent=ext, aspect='auto', cmap='viridis')
    axs[0].set(title='mean (ln)', xlabel='theta_h (deg)', ylabel='theta_d (deg)'); fig.colorbar(im, ax=axs[0], shrink=0.8)
    for j in range(nshow):
        b = _slice_image(model, model.basis[:, j] * np.sqrt(model.eigvals[j] / n), ipd); v = np.nanmax(np.abs(b))
        im = axs[j + 1].imshow(b, origin='lower', extent=ext, aspect='auto', cmap='RdBu_r', vmin=-v, vmax=v)
        axs[j + 1].set(title=f'PC{j + 1} x std ({model.evr[j]:.1%})', xlabel='theta_h (deg)'); fig.colorbar(im, ax=axs[j + 1], shrink=0.8)
    fig.suptitle(f'mean and principal components, channel mean, phi_d = {np.degrees(model.phi_d[ipd]):.0f} deg (theta_h axis is MERL squared spacing)')
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, 'merl_pca_basis.png'), dpi=130); plt.close(fig)

    worst = np.argsort(model.loo_rmse[:, K])[::-1][:3] if model.loo_rmse.shape[1] > K else np.argsort(model.recon_rmse[:, K])[::-1][:3]
    best = [int(np.argsort(model.loo_rmse[:, K])[0])] if model.loo_rmse.shape[1] > K else [int(np.argsort(model.recon_rmse[:, K])[0])]
    rows = [int(i) for i in best] + [int(i) for i in worst]
    fig, axs = plt.subplots(len(rows), 3, figsize=(11, 3.1 * len(rows)))
    for r, i in enumerate(rows):
        x = model.mean + Xc[i]; rec = model.mean + model.basis[:, :K] @ model.coef[i, :K]
        A = np.log10(np.maximum(model.unmap(x), 1e-6)).ravel(); B = np.log10(np.maximum(model.unmap(rec), 1e-6)).ravel()
        imgs = [_slice_image(model, A, ipd), _slice_image(model, B, ipd)]
        vmin, vmax = np.nanmin(imgs[0]), np.nanmax(imgs[0])
        for c, (img, ttl) in enumerate(zip(imgs, (f'{model.names[i]}  log10 f', f'reconstruction k={K}'))):
            im = axs[r, c].imshow(img, origin='lower', extent=ext, aspect='auto', cmap='magma', vmin=vmin, vmax=vmax); axs[r, c].set(title=ttl); fig.colorbar(im, ax=axs[r, c], shrink=0.8)
        d = imgs[1] - imgs[0]; v = min(np.nanmax(np.abs(d)), 1.0)
        im = axs[r, 2].imshow(d, origin='lower', extent=ext, aspect='auto', cmap='RdBu_r', vmin=-v, vmax=v)
        axs[r, 2].set(title=f'log10 difference  (LOO RMSE {model.loo_rmse[i, K]:.3f})' if model.loo_rmse.shape[1] > K else 'log10 difference'); fig.colorbar(im, ax=axs[r, 2], shrink=0.8)
    for a in axs[-1]: a.set_xlabel('theta_h (deg)')
    for a in axs[:, 0]: a.set_ylabel('theta_d (deg)')
    fig.suptitle(f'best and worst materials (by leave-one-out RMSE), phi_d = {np.degrees(model.phi_d[ipd]):.0f} deg'); fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'merl_pca_recon.png'), dpi=130); plt.close(fig)
