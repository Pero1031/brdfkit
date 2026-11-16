"""
Streaming two-part PCA of a BRDFTable (one phi_d slice at a time, so a 13 GB hpBRDF never has to fit in memory).

Per angular bin:
  part A  x1 = log(max(m00, eps))                      (n_ch,)      -> k1 coefficients
  part B  x2 = clip(m_ij / m00, -1, 1), ij != 00       (n_ch * 15,) -> k2 coefficients   (polarized tables only)
Decode:   m00 = exp(mean1 + B1 c1);  m_ij = m00 * (mean2 + B2 c2)
Unmeasured bins (NaN) are excluded from the fit and, in the coefficient table, filled with the nearest measured bin of the
same phi_d slice; the `valid` mask is stored alongside.

Output .npz keys: coef_logm00 [P,D,H,k1], coef_ratio [P,D,H,k2], valid [P,D,H], mean_logm00, basis_logm00 [C,k1], evr_logm00,
mean_ratio, basis_ratio [C*15,k2], evr_ratio, ratio_index, phi_d, theta_d, theta_h, channels, channel_kind, polarized, meta (json).
A sibling *_basis_full.npz keeps every eigenvector for k sweeps.
"""
import json
import os
import time
import numpy as np
from scipy import ndimage
from ..table import BRDFTable

EPS = 1e-4


def ratio_index(n_ch):
    """indices into flatten(n_ch, 16) of all Mueller elements except m00"""
    return np.arange(n_ch * 16).reshape(n_ch, 16)[:, 1:].ravel()


def features(sl, polarized, eps=EPS):
    """slice -> valid mask [D,H], x1 (n, n_ch), x2 (n, n_ch*15 or 0)"""
    nd = 3 if polarized else 1
    flat = sl.reshape(sl.shape[:-nd] + (-1,))
    valid = np.isfinite(flat).all(-1)
    V = flat[valid].astype(np.float64)
    if polarized:
        n_ch = sl.shape[-3]
        m00 = np.maximum(V[:, 0::16], eps)
        x2 = np.clip(V[:, ratio_index(n_ch)] / np.repeat(m00, 15, axis=1), -1.0, 1.0)
    else:
        m00 = np.maximum(V, eps)
        x2 = np.zeros((len(V), 0))
    return valid, np.log(m00), x2


class CovAcc:
    def __init__(self, d):
        self.d = d; self.n = 0; self.s = np.zeros(d); self.S = np.zeros((d, d))

    def add(self, X):
        if len(X) == 0 or self.d == 0:
            return
        self.n += len(X); self.s += X.sum(0); self.S += X.T @ X

    def pca(self):
        if self.d == 0:
            return np.zeros(0), np.zeros((0, 0)), np.zeros(0)
        mean = self.s / self.n
        cov = self.S / self.n - np.outer(mean, mean)
        w, V = np.linalg.eigh(cov)
        order = np.argsort(w)[::-1]
        w = np.maximum(w[order], 0); V = V[:, order]
        return mean, V, w / max(w.sum(), 1e-300)


def compress(table: BRDFTable, out, k1=8, k2=32, dtype='float16', step=1, log=print):
    """Fit the two PCAs over the table and write the coefficient table to `out` (.npz). Returns the npz path."""
    n_phi, n_thd, n_thh, n_ch = table.n_phi, table.n_thd, table.n_thh, table.n_ch
    if not table.polarized:
        k2 = 0
    k1 = min(k1, n_ch); k2 = min(k2, n_ch * 15)
    log(table.describe())

    t0 = time.time(); A1 = CovAcc(n_ch); A2 = CovAcc(n_ch * 15 if table.polarized else 0)
    for i, sl in table.slices(step):
        _, x1, x2 = features(sl, table.polarized)
        A1.add(x1); A2.add(x2)
        if i % 40 == 0:
            log(f'pass1 {i}/{n_phi}  {time.time()-t0:.0f}s')
    mean1, V1, evr1 = A1.pca(); mean2, V2, evr2 = A2.pca()
    log('pass1 done %.0fs, samples %d' % (time.time() - t0, A1.n))
    for k in (1, 2, 4, 8, 16, 32, 64):
        line = f'  cum.var  logm00 k={k:3d}: {evr1[:k].sum():.6f}'
        if table.polarized:
            line += f'   ratio k={k:3d}: {evr2[:k].sum():.6f}'
        log(line)
    base = out[:-4] if out.endswith('.npz') else out
    np.savez(base + '_basis_full.npz', mean_logm00=mean1, basis_logm00=V1, evr_logm00=evr1,
             mean_ratio=mean2, basis_ratio=V2, evr_ratio=evr2)

    B1 = V1[:, :k1]; B2 = V2[:, :k2]
    dt = np.dtype(dtype)
    C1 = np.zeros((n_phi, n_thd, n_thh, k1), dt); C2 = np.zeros((n_phi, n_thd, n_thh, k2), dt)
    valid = np.zeros((n_phi, n_thd, n_thh), bool); empty = []
    t0 = time.time()
    for i, sl in table.slices():
        v, x1, x2 = features(sl, table.polarized)
        valid[i] = v
        c1 = np.zeros((n_thd, n_thh, k1)); c2 = np.zeros((n_thd, n_thh, k2))
        if v.any():
            c1[v] = (x1 - mean1) @ B1
            if k2:
                c2[v] = (x2 - mean2) @ B2
            if not v.all():
                idx = ndimage.distance_transform_edt(~v, return_distances=False, return_indices=True)
                c1 = c1[idx[0], idx[1]]; c2 = c2[idx[0], idx[1]]
        else:
            empty.append(i)
        C1[i] = c1; C2[i] = c2
        if i % 40 == 0:
            log(f'pass2 {i}/{n_phi}  {time.time()-t0:.0f}s')
    log('pass2 done %.0fs' % (time.time() - t0))

    meta = dict(table.meta, k1=k1, k2=k2, eps=EPS, coef_dtype=str(dtype), empty_slices=empty,
                layout='coef_*[phi_d, theta_d, theta_h, k]; Mueller flatten index = ch*16 + row*4 + col',
                decode='m00 = exp(mean_logm00 + basis_logm00 @ c1); m_ij = m00 * (mean_ratio + basis_ratio @ c2)[ratio_index]')
    meta = {k: (v if isinstance(v, (str, int, float, bool, list, dict, type(None))) else str(v)) for k, v in meta.items()}
    np.savez(out, coef_logm00=C1, coef_ratio=C2, valid=valid,
             mean_logm00=mean1.astype(np.float32), basis_logm00=B1.astype(np.float32), evr_logm00=evr1[:k1].astype(np.float32),
             mean_ratio=mean2.astype(np.float32), basis_ratio=B2.astype(np.float32), evr_ratio=evr2[:k2].astype(np.float32),
             ratio_index=ratio_index(n_ch).astype(np.int32),
             phi_d=table.phi_d.astype(np.float32), theta_d=table.theta_d.astype(np.float32), theta_h=table.theta_h.astype(np.float32),
             channels=table.channels.astype(np.float32), channel_kind=table.channel_kind, polarized=table.polarized,
             meta=json.dumps(meta))
    log(f'wrote {out}  {os.path.getsize(out)/2**20:.1f} MB  (coef {(C1.nbytes + C2.nbytes)/2**20:.1f} MB)')
    return out


def evaluate(table: BRDFTable, npz, slices=None, sweep=False, log=print):
    """Reconstruction error of the coefficient table against the source, on measured bins. Returns a dict of metrics."""
    from ..decode import decode
    with np.load(npz) as zf:
        z = dict(zf)
    polarized = bool(z['polarized'])
    if slices is None:
        slices = list(range(0, table.n_phi, max(1, table.n_phi // 12)))
    log('k1 =', z['basis_logm00'].shape[1], ' k2 =', z['basis_ratio'].shape[1], ' coef dtype', z['coef_logm00'].dtype)
    err_m00 = []; err_log = []; err_ratio = []
    for i in slices:
        sl = table.slice(i); v = z['valid'][i]
        R = decode(z['coef_logm00'][i], z['coef_ratio'][i], z)
        o = sl[v]; r = R[v]
        m00o = np.maximum(table.m00(o), EPS); m00r = table.m00(r)
        err_m00.append(np.sqrt(np.mean((m00r - m00o) ** 2)) / np.sqrt(np.mean(m00o ** 2)))
        err_log.append(np.sqrt(np.mean((np.log(m00r) - np.log(m00o)) ** 2)))
        rel = np.abs(m00r - m00o) / m00o
        line = (f'phi_d[{i:3d}]={z["phi_d"][i]:+.2f}  valid {v.mean():.2f}  m00 relRMSE {err_m00[-1]:.4f}  '
                f'logRMSE {err_log[-1]:.4f}  |rel| median {np.median(rel):.4f} p99 {np.percentile(rel, 99):.4f}')
        if polarized:
            ro = np.clip(o / m00o[..., None, None], -1, 1); rr = r / m00r[..., None, None]
            err_ratio.append(np.sqrt(np.mean((rr - ro) ** 2, axis=(0, 1))))
        log(line)
    out = dict(m00_rel_rmse=float(np.mean(err_m00)), log_rmse=float(np.mean(err_log)))
    log('\nmean over slices:  m00 relRMSE %.4f   log m00 RMSE %.4f' % (out['m00_rel_rmse'], out['log_rmse']))
    if polarized:
        out['ratio_rmse'] = np.mean(err_ratio, 0)
        log('RMSE of normalized Mueller elements (m_ij/m00), averaged over slices:')
        log(np.array2string(out['ratio_rmse'], precision=4, suppress_small=True))
    if sweep:
        with np.load(npz[:-4] + '_basis_full.npz') as f_:
            zf = dict(f_)
        i = slices[len(slices) // 2]
        v, x1, x2 = features(table.slice(i), polarized)
        m00o = np.exp(x1)
        log(f'\nk sweep on phi_d[{i}]  (projection of exact features, float64):')
        log('  log m00:  k   relRMSE(m00)   logRMSE')
        for k in (1, 2, 4, 8, 16, 32):
            if k > zf['basis_logm00'].shape[1]: break
            B = zf['basis_logm00'][:, :k]
            rec = np.exp(zf['mean_logm00'] + ((x1 - zf['mean_logm00']) @ B) @ B.T)
            log(f'           {k:3d}   {np.sqrt(np.mean((rec-m00o)**2))/np.sqrt(np.mean(m00o**2)):.5f}      {np.sqrt(np.mean((np.log(rec)-x1)**2)):.5f}')
        if polarized:
            log('  ratios :  k   RMSE(all 15 elems)   max elem RMSE')
            for k in (4, 8, 16, 32, 64, 128):
                if k > zf['basis_ratio'].shape[1]: break
                B = zf['basis_ratio'][:, :k]
                rec = zf['mean_ratio'] + ((x2 - zf['mean_ratio']) @ B) @ B.T
                e = (rec - x2).reshape(len(x2), -1, 15)
                log(f'           {k:3d}   {np.sqrt(np.mean(e**2)):.5f}              {np.sqrt(np.mean(e ** 2, axis=(0, 1))).max():.5f}')
    return out
