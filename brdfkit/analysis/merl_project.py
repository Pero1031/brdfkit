"""
Project a measured table (hpBRDF aluminium, 68 wavelengths, polarized) onto the inter-material PCA basis of the MERL database
(`brdfkit merl-pca`), to see where the material sits among the 100 MERL materials and how well that basis represents it.

Step 1, `merl_grid_rgb`: table -> unpolarized RGB table on the MERL grid (180 x 90 x 90).
  * Grid. The hpBRDF axes (phi_d every 1 deg over [-180, 180], theta_d every 1 deg, theta_h with the same squared spacing as MERL)
    contain every MERL grid point to float precision, so no interpolation is done: bins are picked by index. A table whose
    axes do not contain the MERL grid is rejected.
  * Colour. m00 of the 68 bands -> linear sRGB with the band weights of the renderer (`brdfkit.color`, same as MeasuredHpBRDF.hpp).
  * phi_d. MERL stores phi_d in [0, 180) deg. For an isotropic, reciprocal, mirror-symmetric surface f(phi) = f(-phi) = f(phi + 180)
    = f(180 - phi); sym='full' averages the table slices of these equivalent angles into each MERL bin (and records their relative
    spread, a direct measure of how well the measurement obeys the symmetries); sym='none' takes the slice at the same angle only.
Step 2, `analyse`: least-squares coefficients on the first k components (the table is unmeasured in a thin band of bins at the
horizon that MERL measured, so the projection is masked), residuals for k = 1..K, distances to every MERL material in the
coefficient space, and a comparison with the MERL aluminium when given.
`run_report` does both for one or more saved `MaterialPCA` models and writes reports and figures.
"""
import glob
import json
import os
import time
import numpy as np
from ..table import BRDFTable, RGB, WAVELENGTH
from ..color import rgb_band_weights
from ..io import merl
from .merl_pca import MaterialPCA

TOL = 1e-4                       # rad, axis matching tolerance (the files store float32)


# ----------------------------------------------------------------------------- table -> MERL grid
def _wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def _fold(phi):
    """canonical phi_d in [0, pi/2] under phi -> -phi (mirror) and phi -> phi + pi (reciprocity)"""
    p = np.abs(_wrap(phi))
    return np.minimum(p, np.pi - p)


def _index_of(src, dst, name):
    idx = np.abs(src[:, None] - dst[None, :]).argmin(0)
    err = np.abs(src[idx] - dst).max()
    if err > TOL:
        raise ValueError(f'{name} axis of the table does not contain the MERL grid (max mismatch {err:.1e} rad); resampling is not implemented')
    return idx


def merl_grid_rgb(table, sym='full', log=print):
    """-> (BRDFTable RGB unpolarized on the MERL grid, info) with info = dict(spread [180, 90, 90], n_slices, seconds)"""
    if sym not in ('full', 'none'):
        raise ValueError("sym must be 'full' or 'none'")
    phi_m, thd_m, thh_m = merl.grid_axes()
    itd, ith = _index_of(table.theta_d, thd_m, 'theta_d'), _index_of(table.theta_h, thh_m, 'theta_h')
    if table.channel_kind == WAVELENGTH:
        weights = rgb_band_weights(table.channels)
    elif table.n_ch == 3:
        weights = None
    else:
        raise ValueError('expected a wavelength table or an RGB table')
    pj, pb = table.phi_d[:, None], phi_m[None, :]
    A = (np.abs(_fold(pj) - _fold(pb)) if sym == 'full' else np.abs(_wrap(pj - pb))) < TOL        # [slice, MERL bin]
    if table.n_phi > 1 and abs(abs(table.phi_d[-1] - table.phi_d[0]) - 2 * np.pi) < TOL:
        A[-1] = False                                                                 # -180 and +180 deg are the same direction
    if not A.any(0).all():
        raise ValueError(f'phi_d of the table does not cover MERL bins {np.nonzero(~A.any(0))[0][:5]}...')
    shape = (len(phi_m), len(thd_m), len(thh_m), 3)
    S, S2, N = np.zeros(shape), np.zeros(shape), np.zeros(shape[:3])
    rows = np.nonzero(A.any(1))[0]
    t0 = time.time()
    for c, j in enumerate(rows):
        m00 = table.m00(table.slice(j))[np.ix_(itd, ith)].astype(np.float64)         # [90, 90, n_ch]
        rgb = m00 if weights is None else m00 @ weights
        fin = np.isfinite(rgb).all(-1)
        rgb = np.where(fin[..., None], rgb, 0.0)
        for b in np.nonzero(A[j])[0]:
            S[b] += rgb; S2[b] += rgb * rgb; N[b] += fin
        if c % 60 == 0 or c == len(rows) - 1:
            log(f'  slice {c + 1}/{len(rows)}  {time.time() - t0:.0f}s')
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = S / N[..., None]
        std = np.sqrt(np.maximum(S2 / N[..., None] - mean ** 2, 0))
        spread = np.where(N > 1, (std / np.maximum(mean, 1e-12)).mean(-1), np.nan)
    values = np.maximum(mean, 0.0).astype(np.float32)                                 # the renderer clamps the baked RGB at 0
    meta = dict(source=str(table.meta.get('source', '')), sym=sym, kind='m00 -> linear sRGB on the MERL grid',
                slices_used=int(len(rows)), slices_per_bin=[int(N[N > 0].min()), int(N.max())])
    out = BRDFTable(phi_m, thd_m, thh_m, np.arange(3.0), RGB, False, meta, values=values)
    return out, dict(spread=spread.astype(np.float32), n_slices=int(len(rows)), seconds=time.time() - t0)


def save_grid(path, table, info):
    np.savez(path, values=table.values, spread=info['spread'], meta=json.dumps(table.meta))
    return path


def load_grid(path):
    with np.load(path, allow_pickle=False) as z:
        meta = json.loads(str(z['meta']))
        v, spread = z['values'], z['spread']
    phi_m, thd_m, thh_m = merl.grid_axes()
    return BRDFTable(phi_m, thd_m, thh_m, np.arange(3.0), RGB, False, meta, values=v), dict(spread=spread)


def cached_grid(source, cache, sym='full', recompute=False, log=print):
    """merl_grid_rgb of the table file `source`, stored in / read from `cache` (the 13 GB read takes a minute or two)"""
    from ..io import load
    if not recompute and os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(source):
        t, info = load_grid(cache)
        if t.meta.get('sym') == sym:
            log(f'using cached MERL-grid table {cache}')
            return t, info
    t, info = merl_grid_rgb(load(source), sym, log)
    save_grid(cache, t, info)
    return t, info


# ----------------------------------------------------------------------------- projection
def linear_rel_rmse(model, values, coef, k):
    """relative RMSE of the cos-weighted linear BRDF f cos_i cos_o between `values` and its k-component reconstruction,
    over the bins measured in the table (same quantity as `linear_error` in merl_pca)"""
    sp, sd, sh = model.stride
    v = np.asarray(values)[::sp, ::sd, ::sh]
    rec = model.reconstruct(coef, k)
    ok = model.valid & np.isfinite(v).all(-1) & np.isfinite(rec).all(-1)
    w = model.cosw[ok][:, None]
    a, b = v[ok] * w, rec[ok] * w
    return float(np.sqrt(np.mean((b - a) ** 2)) / np.sqrt(np.mean(a ** 2)))


def analyse(model, table, k=16, reference=None, ref_name=None, ks=(1, 2, 3, 4, 5, 6, 8, 10, 12, 16)):
    """project `table` (RGB on the MERL grid) on `model`; `reference` = optional MERL table of a comparison material"""
    k = min(k, model.k)
    ks = [kk for kk in ks if kk <= k]
    coef, info = model.project_masked(table.values, k)
    res = dict(k=k, coef=coef, info=info, ks=ks)
    per_k = [model.project_masked(table.values, kk) for kk in ks]
    res['rmse_k'] = [p[1]['rmse'] for p in per_k]
    res['lin_k'] = [linear_rel_rmse(model, table.values, p[0], kk) for kk, p in zip(ks, per_k)]
    C = model.coef[:, :k]
    res['dist'] = np.linalg.norm(C - coef, axis=1)                                   # [n], in feature-space norm
    res['dist8'] = np.linalg.norm(C[:, :8] - coef[:8], axis=1)
    res['z'] = coef / C.std(0)
    if model.loo_rmse.shape[1] > k:
        res['loo'] = model.loo_rmse[:, k]
    if reference is not None:
        name = ref_name if ref_name in model.names else None
        rc = model.coef[model.names.index(name), :k] if name else model.project_masked(reference.values, k)[0]
        xh, xr = model.features(table.values), model.features(reference.values)
        both = np.isfinite(xh) & np.isfinite(xr)
        vh, vr = table.values.reshape(-1, 3), reference.values.reshape(-1, 3)
        ok = np.isfinite(vh).all(-1) & np.isfinite(vr).all(-1) & (vh > 1e-3).all(-1) & (vr > 1e-3).all(-1)
        pca_dist = float(np.linalg.norm(rc - coef))
        res['ref'] = dict(name=name or ref_name, coef=rc, in_fit=name is not None,
                          raw_rmse=float(np.sqrt(np.mean((xh[both] - xr[both]) ** 2))), pca_dist=pca_dist,
                          log_ratio_rgb=np.median(np.log(vh[ok] / vr[ok]), 0), rank=int((res['dist'] < pca_dist).sum()) + 1)
        i = model.names.index(name) if name else None
        res['ref']['rmse_fit'] = float(model.recon_rmse[i, k]) if name else None
        res['ref']['rmse_loo'] = float(model.loo_rmse[i, k]) if name and model.loo_rmse.shape[1] > k else None
    return res


# ----------------------------------------------------------------------------- report
def _report(model, res, grid_info, source, say):
    k, n, D = res['k'], model.n, model.dim
    info = res['info']
    say(f'=== {source} on the MERL PCA ({model.mapping} mapping, eps {model.eps:g}, stride {model.stride}, {n} materials, D = {D}) ===')
    sp = grid_info.get('spread')
    if sp is not None and np.isfinite(sp).any():
        s = sp[np.isfinite(sp)]
        say(f'phi_d symmetry: relative spread of the equivalent slices (phi, -phi, phi+180, 180-phi) per bin: median {np.median(s):.3f}, '
            f'p90 {np.percentile(s, 90):.3f}')
    say(f'measured features inside the model mask: {info["n_used"]} of {info["n_total"]} ({info["n_used"] / info["n_total"]:.1%});'
        f' condition number of the restricted basis {info["cond"]:.2f}')
    say('\nk   feature RMSE (ln)   linear rel. RMSE of f cos cos')
    for kk, r, l in zip(res['ks'], res['rmse_k'], res['lin_k']):
        say(f'{kk:3d}   {r:8.3f}            {l:8.3f}')
    if 'loo' in res:
        l = res['loo']
        say(f'MERL materials held out of the fit, k = {k}: median {np.median(l):.3f}, p90 {np.percentile(l, 90):.3f}, max {l.max():.3f}  (feature RMSE, same unit)')
    say(f'\ncoefficients (k = {k}) and their position in the MERL distribution (z = coefficient / std over the {n} materials):')
    say('comp  coef      z     ' + ('MERL ' + res['ref']['name'] + ' coef   z' if 'ref' in res else ''))
    sd = model.coef[:, :k].std(0)
    for j in range(k):
        line = f'{j + 1:4d}  {res["coef"][j]:8.1f}  {res["z"][j]:5.2f}'
        if 'ref' in res:
            line += f'   {res["ref"]["coef"][j]:8.1f}        {res["ref"]["coef"][j] / sd[j]:5.2f}'
        say(line)
    say(f'\nnearest MERL materials in the {k}-d coefficient space (distance in feature-norm units; RMS per feature = distance / sqrt(D) = distance / {np.sqrt(D):.0f}):')
    order = np.argsort(res['dist'])
    for i in order[:8]:
        say(f'  {model.names[i]:26s} {res["dist"][i]:8.1f}  (RMS {res["dist"][i] / np.sqrt(D):.3f})   first 8 comps only: {res["dist8"][i]:8.1f}')
    say('  farthest: ' + ', '.join(f'{model.names[i]} {res["dist"][i]:.0f}' for i in order[::-1][:3]))
    if 'ref' in res:
        r = res['ref']
        say(f'\ncomparison with MERL "{r["name"]}" ({"in the fit" if r["in_fit"] else "not in the fit"}):')
        say(f'  raw feature RMSE hpBRDF vs MERL, all measured bins (no PCA): {r["raw_rmse"]:.3f}')
        say(f'  distance in the {k}-d coefficient space: {r["pca_dist"]:.1f} (RMS {r["pca_dist"] / np.sqrt(D):.3f}); '
            f'rank among the MERL materials by distance to the hpBRDF material: {r["rank"]}')
        say('  median ln(hpBRDF / MERL) per channel R, G, B: ' + ', '.join(f'{v:+.3f}' for v in r['log_ratio_rgb']) +
            '   (' + ', '.join(f'x{np.exp(v):.2f}' for v in r['log_ratio_rgb']) + ')')
        if r['rmse_fit'] is not None:
            say(f'  the MERL material itself at k = {k}: reconstruction RMSE {r["rmse_fit"]:.3f}' +
                (f', leave-one-out RMSE {r["rmse_loo"]:.3f}' if r['rmse_loo'] is not None else ''))


# ----------------------------------------------------------------------------- figures
def _log10_image(values, model, ipd):
    """linear values on the unstrided grid [P, T, H, 3] or the strided grid -> channel-mean log10 image [T', H'] at phi_d slice ipd"""
    v = np.asarray(values)
    if v.shape[:3] != model.grid_shape:
        sp, sd, sh = model.stride
        v = v[::sp, ::sd, ::sh]
    img = np.log10(np.maximum(v[ipd], 1e-6)).mean(-1)
    img[~model.valid[ipd]] = np.nan
    return img


def _plots(model, res, table, reference, source, out_dir, tag):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    k, n = res['k'], model.n
    C = model.coef; c = res['coef']
    order = np.argsort(res['dist'])
    ref = res.get('ref')
    ri = model.names.index(ref['name']) if ref and ref['in_fit'] else None

    # --- scatter
    fig, axs = plt.subplots(1, 2, figsize=(17, 8))
    for a, (i, j) in zip(axs, ((0, 1), (2, 3))):
        a.scatter(C[:, i], C[:, j], c='0.65', s=22, edgecolor='k', lw=0.3, label='MERL materials')
        for m in set(order[:6].tolist()) | ({ri} if ri is not None else set()):
            a.annotate(model.names[m], (C[m, i], C[m, j]), fontsize=7, xytext=(3, 3), textcoords='offset points')
        if ri is not None:
            a.scatter(C[ri, i], C[ri, j], marker='D', s=70, c='tab:blue', edgecolor='k', label=f'MERL {ref["name"]}')
            a.plot([C[ri, i], c[i]], [C[ri, j], c[j]], 'r-', lw=0.8)
        a.scatter(c[i], c[j], marker='*', s=260, c='tab:red', edgecolor='k', label=source, zorder=5)
        a.set(xlabel=f'PC{i + 1} ({model.evr[i]:.1%})', ylabel=f'PC{j + 1} ({model.evr[j]:.1%})'); a.grid(alpha=0.3)
        a.axhline(0, color='gray', lw=0.5); a.axvline(0, color='gray', lw=0.5)
    axs[0].legend(fontsize=9)
    fig.suptitle(f'{source} in the MERL principal-component space ({model.mapping} mapping, k = {k})')
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, f'{tag}_scatter.png'), dpi=130); plt.close(fig)

    # --- slices at phi_d = 90 deg, difference maps and a theta_h profile
    ipd = model.grid_shape[0] // 2
    ext = [np.degrees(model.theta_h[0]), np.degrees(model.theta_h[-1]), np.degrees(model.theta_d[0]), np.degrees(model.theta_d[-1])]
    rec = model.reconstruct(c, k)
    imgs = [(_log10_image(table.values, model, ipd), f'{source}  log10 f'), (_log10_image(rec, model, ipd), f'projection, k = {k}')]
    if reference is not None:
        imgs.append((_log10_image(reference.values, model, ipd), f'MERL {ref["name"]}  log10 f'))
    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    lo, hi = np.nanmin(imgs[0][0]), np.nanmax(imgs[0][0])
    for a, (img, ttl) in zip(axs[0], imgs):
        im = a.imshow(img, origin='lower', extent=ext, aspect='auto', cmap='magma', vmin=lo, vmax=hi)
        a.set(title=ttl, xlabel='theta_h (deg, squared MERL spacing)', ylabel='theta_d (deg)'); fig.colorbar(im, ax=a, shrink=0.8)
    for a in axs[0][len(imgs):]:
        a.axis('off')
    diffs = [(imgs[1][0] - imgs[0][0], f'projection - {source}')]
    if reference is not None:
        diffs.append((imgs[2][0] - imgs[0][0], f'MERL {ref["name"]} - {source}'))
    for a, (d, ttl) in zip(axs[1], diffs):
        v = min(np.nanmax(np.abs(d)), 1.0)
        im = a.imshow(d, origin='lower', extent=ext, aspect='auto', cmap='RdBu_r', vmin=-v, vmax=v)
        a.set(title=ttl + ' (log10)', xlabel='theta_h (deg)', ylabel='theta_d (deg)'); fig.colorbar(im, ax=a, shrink=0.8)
    itd = int(np.argmin(np.abs(np.degrees(model.theta_d) - 30)))
    a = axs[1][2]; th = np.degrees(model.theta_h)
    sp, sd, sh = model.stride
    ok = model.valid[ipd, itd]
    prof = lambda v: np.where(ok[:, None], np.asarray(v)[::sp, ::sd, ::sh][ipd, itd], np.nan)
    for ch, col in enumerate(('tab:red', 'tab:green', 'tab:blue')):
        a.semilogy(th, prof(table.values)[:, ch], '-', c=col, lw=2, label='hpBRDF' if ch == 0 else None)
        a.semilogy(th, np.where(ok, rec[ipd, itd][:, ch], np.nan), '--', c=col, lw=1.2, label='projection' if ch == 0 else None)
        if reference is not None:
            a.semilogy(th, prof(reference.values)[:, ch], ':', c=col, lw=1.2, label=('MERL ' + ref['name']) if ch == 0 else None)
    a.set(title=f'profile at theta_d = {np.degrees(model.theta_d[itd]):.0f} deg, R G B', xlabel='theta_h (deg)', ylabel='f (1/sr)')
    a.grid(alpha=0.3, which='both'); a.legend(fontsize=8)
    fig.suptitle(f'phi_d = {np.degrees(model.phi_d[ipd]):.0f} deg, channel mean ({model.mapping} mapping)'); fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f'{tag}_slices.png'), dpi=130); plt.close(fig)

    # --- residual versus k, coefficient z-scores
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.6))
    a = axs[0]
    if model.loo_rmse.shape[1] > 1:
        kk = np.arange(1, min(n - 1, model.loo_rmse.shape[1]))
        a.fill_between(kk, np.percentile(model.loo_rmse[:, kk], 10, 0), np.percentile(model.loo_rmse[:, kk], 90, 0), color='0.85', label='MERL leave-one-out, p10-p90')
        a.plot(kk, np.median(model.loo_rmse[:, kk], 0), '-', c='0.5', label='MERL leave-one-out, median')
        if ri is not None:
            a.plot(kk, model.loo_rmse[ri, kk], '-', c='tab:blue', label=f'MERL {ref["name"]}, leave-one-out')
            a.plot(kk, model.recon_rmse[ri, kk], '--', c='tab:blue', label=f'MERL {ref["name"]}, in the fit')
    a.plot(res['ks'], res['rmse_k'], 'o-', c='tab:red', lw=2, label=source)
    a.set(xscale='log', xlabel='k', ylabel='feature RMSE (ln units)', title='residual after projection', xlim=(1, max(res['ks']) * 1.1))
    a.grid(alpha=0.3); a.legend(fontsize=8)
    a = axs[1]; x = np.arange(1, k + 1); w = 0.4
    a.bar(x - (w / 2 if ref else 0), res['z'], w, color='tab:red', label=source)
    if ref:
        a.bar(x + w / 2, ref['coef'] / C[:, :k].std(0), w, color='tab:blue', label=f'MERL {ref["name"]}')
    a.set(xlabel='component', ylabel='coefficient / std over MERL', title='position in the MERL distribution (z-score)', xticks=x)
    a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, f'{tag}_error.png'), dpi=130); plt.close(fig)


# ----------------------------------------------------------------------------- driver
def run_report(source, model_paths, out_dir, k=16, sym='full', reference=None, recompute=False, plot=True, log=print):
    """source: hpBRDF-like table file; model_paths: saved MaterialPCA .npz files (one report each); reference: MERL .binary of the
    comparison material (e.g. data/merl/aluminium.binary). Writes <out_dir>/<mapping>/{report.txt, projection.json, *.png}
    and the MERL-grid RGB table <out_dir>/<name>_merl_grid_<sym>.npz (reused on the next run)."""
    from ..io import load
    os.makedirs(out_dir, exist_ok=True)
    name = os.path.splitext(os.path.basename(source))[0]
    table, ginfo = cached_grid(source, os.path.join(out_dir, f'{name}_merl_grid_{sym}.npz'), sym, recompute, log)
    ref = load(reference) if reference else None
    ref_name = os.path.splitext(os.path.basename(reference))[0] if reference else None
    summary = {}
    for mp in model_paths:
        model = MaterialPCA.load(mp)
        sub = os.path.join(out_dir, model.mapping); os.makedirs(sub, exist_ok=True)
        res = analyse(model, table, k, ref, ref_name)
        lines = []

        def say(*a):
            s = ' '.join(str(x) for x in a); log(s); lines.append(s)
        say(f'model: {mp}')
        _report(model, res, ginfo, name, say)
        open(os.path.join(sub, 'report.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
        js = dict(source=source, model=mp, mapping=model.mapping, k=res['k'], sym=sym, coef=res['coef'].tolist(), info=res['info'],
                  ks=res['ks'], rmse_k=res['rmse_k'], linear_rel_rmse_k=res['lin_k'],
                  nearest=[dict(name=model.names[i], dist=float(res['dist'][i])) for i in np.argsort(res['dist'])[:10]])
        if 'ref' in res:
            r = res['ref']
            js['reference'] = dict(name=r['name'], in_fit=r['in_fit'], coef=r['coef'].tolist(), raw_rmse=r['raw_rmse'], pca_dist=r['pca_dist'],
                                   rank=r['rank'], log_ratio_rgb=r['log_ratio_rgb'].tolist())
        json.dump(js, open(os.path.join(sub, 'projection.json'), 'w'), indent=1)
        if plot:
            _plots(model, res, table, ref, name, sub, 'hpbrdf_on_merl')
        summary[model.mapping] = res
    return summary


# ----------------------------------------------------------------------------- comparison of several runs
def grid_stats(values, iphi=90, itd=30):
    """peak shape and floor of a table on the MERL grid (green channel): peak_max over the grid, peak_core = f at theta_h = 0 for
    phi_d = iphi deg and theta_d = itd deg, flat_ratio = f(theta_h = 0.28 deg) / f(0) there (1 = flat top), floor = median of f over
    theta_h >= 40 deg and theta_d < 80 deg (far from the specular lobe)"""
    g = values[iphi, itd, :, 1]
    return dict(peak_max=float(np.nanmax(values[..., 1])), peak_core=float(g[0]), flat_ratio=float(g[5] / g[0]),
                floor=float(np.nanmedian(values[:, :80, 60:, 1])))


def compare_runs(run_dirs, say=print):
    """summarize several `merl-project` output folders (one table each): peak, floor, phi_d spread, and per mapping the residual,
    the most extreme component, the nearest MERL materials and the comparison with the reference material"""
    for d in run_dirs:
        npz = glob.glob(os.path.join(d, '*_merl_grid_full.npz'))
        if not npz:
            raise ValueError(f'{d}: no *_merl_grid_full.npz (run merl-project with --sym full first)')
        name = os.path.basename(npz[0]).split('_merl_grid')[0]
        t, info = load_grid(npz[0])
        st, sp = grid_stats(t.values), info['spread'][np.isfinite(info['spread'])]
        say(f'\n== {name}: max {st["peak_max"]:.2f}, f(theta_h = 0) {st["peak_core"]:.2f} at phi_d 90, theta_d 30 deg, flat-top ratio {st["flat_ratio"]:.3f}, '
            f'floor {st["floor"]:.4f}, phi_d spread median {np.median(sp):.3f}')
        for mapping in ('log', 'nielsen'):
            pj = json.load(open(os.path.join(d, mapping, 'projection.json')))
            model = MaterialPCA.load(pj['model'])
            c = np.array(pj['coef']); z = c / model.coef[:, :len(c)].std(0)
            ref = pj.get('reference', {})
            say(f'  {mapping:8s} rmse k=8 {pj["rmse_k"][pj["ks"].index(8)]:.3f} k={pj["k"]} {pj["rmse_k"][-1]:.3f}; max |z| {np.abs(z).max():.1f} at PC{np.abs(z).argmax() + 1}, '
                f'z7 {z[6]:+.1f}; nearest {", ".join(x["name"] for x in pj["nearest"][:4])}' +
                (f'; {ref["name"]}: rank {ref["rank"]}, raw RMSE {ref["raw_rmse"]:.2f}, median ratio RGB ' +
                 '/'.join(f'{np.exp(v):.2f}' for v in ref['log_ratio_rgb']) if ref else ''))
