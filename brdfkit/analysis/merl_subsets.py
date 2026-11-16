"""
Does a metal-only PCA basis represent a measured metal better than the basis of all 100 MERL materials?

`brdfkit merl-project` showed that the hpBRDF aluminium sits far from the MERL aluminium in the 100-material PCA and leans on one
component (PC7, the specular core). Two explanations: the all-material basis lacks the variety of metal specular shapes, or the table itself
is unlike any MERL metal. This module separates them by refitting the PCA on subsets of the MERL materials and projecting the table on each:

  * 'metals'           11 bare metals / steels (METALS)
  * 'metals+finishes'  those plus 12 metallic paints and two-layer coatings (FINISHES)
  * each with and without the MERL aluminium (the same material class as the table; without it the basis never saw aluminium)
  * random subsets of the same sizes drawn from the other materials: the control. A metal basis helps only if it beats these.

All fits use the log mapping on the stride-2 grid. A basis of fewer materials has a (slightly) larger common mask, so every residual is
measured on the same bins, those measured in all 100 materials and in the table, and only there: ln RMSE of the log BRDF and the
relative RMSE of the cos-weighted linear BRDF. The coefficients are the masked least-squares fit of `MaterialPCA.project_masked`.
"""
import json
import os
import time
import numpy as np
from ..table import BRDFTable
from ..io import load
from .merl_pca import MaterialPCA, fit, list_materials
from .merl_project import load_grid

EPS = 1e-5
METALS = ('aluminium', 'alum-bronze', 'brass', 'chrome', 'chrome-steel', 'nickel', 'ss440', 'steel', 'tungsten-carbide',
          'black-oxidized-steel', 'grease-covered-steel')
FINISHES = ('gold-metallic-paint', 'gold-metallic-paint2', 'gold-metallic-paint3', 'silver-metallic-paint', 'silver-metallic-paint2',
            'blue-metallic-paint', 'blue-metallic-paint2', 'green-metallic-paint', 'green-metallic-paint2', 'red-metallic-paint',
            'two-layer-gold', 'two-layer-silver')
KS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16)
KMAX = 16


def _strided(t, s=2):
    return BRDFTable(t.phi_d[::s], t.theta_d[::s], t.theta_h[::s], t.channels, t.channel_kind, t.polarized, t.meta,
                     values=np.ascontiguousarray(t.values[::s, ::s, ::s]))


def residuals(model, hp, common, ks):
    """-> (ln RMSE [len(ks)], relative linear RMSE [len(ks)], max |z| of the coefficients at the largest k, its component)
    hp [P, T, H, 3] on the model grid; common [P, T, H] bool"""
    ln, lin = [], []
    for k in ks:
        c, _ = model.project_masked(hp, k)
        rec = model.reconstruct(c, k)
        ok = common & np.isfinite(rec).all(-1)
        a, b = hp[ok].astype(np.float64), rec[ok].astype(np.float64)
        ln.append(float(np.sqrt(np.mean((np.log(np.maximum(a, EPS)) - np.log(np.maximum(b, EPS))) ** 2))))
        w = model.cosw[ok][:, None]
        lin.append(float(np.sqrt(np.mean(((b - a) * w) ** 2)) / np.sqrt(np.mean((a * w) ** 2))))
    k = ks[-1]
    z = model.project_masked(hp, k)[0] / model.coef[:, :k].std(0)
    return np.array(ln), np.array(lin), float(np.abs(z).max()), int(np.abs(z).argmax()) + 1


def run_study(hp_grid, merl_dir, out_dir, all_model=None, n_random=30, seed=0, ablate=False, log=print):
    """hp_grid: the MERL-grid table written by `merl-project` (<out>/<name>_merl_grid_full.npz); merl_dir: the .binary files;
    all_model: saved 100-material merl_pca.npz (log mapping), for the leave-one-out numbers of the all-material basis;
    ablate: also refit the 100-material basis 100 times, each time without one material, and report how the residual of the table changes"""
    os.makedirs(out_dir, exist_ok=True)
    quiet = lambda *a: None
    hp_t, _ = load_grid(hp_grid)
    hp = np.ascontiguousarray(hp_t.values[::2, ::2, ::2])
    t0 = time.time()
    tables = {}
    for i, (name, path) in enumerate(list_materials(merl_dir)):
        tables[name] = _strided(load(path))
    names = sorted(tables)
    log(f'loaded {len(names)} materials in {time.time() - t0:.0f}s')
    for s in (METALS, FINISHES):
        missing = [m for m in s if m not in tables]
        if missing:
            raise ValueError(f'unknown MERL materials: {missing}')

    def fit_set(sel, loo=False):
        return fit([(n, tables[n]) for n in sel], stride=1, mapping='log', k=KMAX, loo=loo, log=quiet)[0]

    base = fit_set(names)
    common = base.valid & np.isfinite(hp).all(-1)                                  # measured in all 100 materials and in the table
    log(f'common bins {common.sum()} of {common.size}; all-material fit done {time.time() - t0:.0f}s')
    nm = lambda s: [m for m in s if m != 'aluminium']
    sets = {'all100': names, 'metals': list(METALS), 'metals, no aluminium': nm(METALS),
            'metals+finishes': list(METALS + FINISHES), 'metals+finishes, no aluminium': nm(METALS + FINISHES)}
    res, models = {}, {}
    for label, sel in sets.items():
        m = base if label == 'all100' else fit_set(sel, loo=True)
        ks = [k for k in KS if k <= m.k]
        ln, lin, zmax, zc = residuals(m, hp, common, ks)
        res[label] = dict(n=len(sel), ks=ks, ln=ln.tolist(), lin=lin.tolist(), zmax=zmax, zcomp=zc, evr=m.evr[:5].tolist())
        if m.loo_rmse.shape[1] > 1:
            res[label]['loo_median_at_k'] = {int(k): float(np.median(m.loo_rmse[:, k])) for k in ks if k < m.loo_rmse.shape[1]}
        models[label] = m
        log(f'{label}: n = {len(sel)}, ln RMSE at k = {ks[-1]}: {ln[-1]:.3f}  ({time.time() - t0:.0f}s)')
    for label in ('metals', 'metals+finishes'):
        models[label].save(os.path.join(out_dir, f'merl_pca_{label.replace("+", "_")}.npz'))
    # leave-one-out of the metals under the all-material basis, same features and mask (saved model), for matched k
    if all_model:
        am = MaterialPCA.load(all_model)
        idx = [am.names.index(n) for n in METALS]
        res['all100']['metals_loo_at_k'] = {int(k): float(np.median(am.loo_rmse[idx, k])) for k in KS if k < am.loo_rmse.shape[1]}
        mm = models['metals']
        res['metals']['metals_loo_at_k'] = {int(k): float(np.median(mm.loo_rmse[:, k])) for k in KS if k < mm.loo_rmse.shape[1]}
    # control: random subsets of the same sizes (never containing the MERL aluminium)
    rng = np.random.default_rng(seed)
    pool = nm(names)
    ctrl = {}
    for size in (len(METALS), len(METALS + FINISHES)):
        runs = []
        for r in range(n_random):
            sel = [str(x) for x in rng.choice(pool, size, replace=False)]
            m = fit_set(sel)
            ks = [k for k in KS if k <= m.k]
            runs.append(residuals(m, hp, common, ks)[:2])
        ctrl[size] = dict(ks=ks, ln=np.array([r[0] for r in runs]), lin=np.array([r[1] for r in runs]))
        log(f'random subsets of {size}: {n_random} fits done ({time.time() - t0:.0f}s)')
    abl = None
    if ablate:
        abl = {}
        for i, n in enumerate(names):
            m = fit_set([x for x in names if x != n])
            abl[n] = residuals(m, hp, common, [8, 16])[0].tolist()
            if i % 20 == 19:
                log(f'ablation {i + 1}/{len(names)}  ({time.time() - t0:.0f}s)')
    _report(res, ctrl, n_random, out_dir, log, abl)
    _plot(res, ctrl, out_dir)
    json.dump(dict(sets=res, random={str(k): dict(ks=v['ks'], ln_mean=v['ln'].mean(0).tolist(), ln_std=v['ln'].std(0).tolist(),
                                                  ln_min=v['ln'].min(0).tolist(), lin_mean=v['lin'].mean(0).tolist()) for k, v in ctrl.items()},
                   metals=list(METALS), finishes=list(FINISHES), n_random=n_random, seed=seed, ablation=abl),
              open(os.path.join(out_dir, 'results.json'), 'w'), indent=1)
    return res, ctrl


def _report(res, ctrl, n_random, out_dir, log, abl=None):
    lines = []

    def say(s=''):
        log(s); lines.append(s)
    say('hpBRDF-derived table projected on PCA bases fitted to subsets of the MERL materials (log mapping, stride 2).')
    say('All residuals on the same bins: measured in all 100 materials and in the table.')
    for key, title in (('ln', 'ln RMSE of the log BRDF'), ('lin', 'relative RMSE of f cos_i cos_o')):
        say(f'\n{title}, by number of components k')
        say(f'{"basis":32s} n   ' + ' '.join(f'{k:>6d}' for k in KS))
        for label, r in res.items():
            vals = dict(zip(r['ks'], r[key]))
            say(f'{label:32s} {r["n"]:3d} ' + ' '.join(f'{vals[k]:6.3f}' if k in vals else '     -' for k in KS))
        for size, c in ctrl.items():
            vals = dict(zip(c['ks'], c[key].mean(0))); sd = dict(zip(c['ks'], c[key].std(0)))
            say(f'{"random " + str(size) + " (mean)":32s} {size:3d} ' + ' '.join(f'{vals[k]:6.3f}' if k in vals else '     -' for k in KS))
            say(f'{"random " + str(size) + " (std)":32s}     ' + ' '.join(f'{sd[k]:6.3f}' if k in sd else '     -' for k in KS))
            if key == 'ln':
                best = dict(zip(c['ks'], c[key].min(0)))
                say(f'{"random " + str(size) + " (best of " + str(n_random) + ")":32s}     ' + ' '.join(f'{best[k]:6.3f}' if k in best else '     -' for k in KS))
    say('\nlargest |z| of the coefficients at the largest k (z = coefficient / std of that component over the basis materials):')
    for label, r in res.items():
        say(f'  {label:32s} k = {r["ks"][-1]:2d}: |z| = {r["zmax"]:5.2f} at component {r["zcomp"]}, first components explain ' +
            ' '.join(f'{v:.3f}' for v in r['evr'][:3]))
    say('\nleave-one-out RMSE of the 11 MERL metals (median over them), same units as the MERL PCA report:')
    for label in ('all100', 'metals'):
        d = res[label].get('metals_loo_at_k')
        if d:
            say(f'  basis {label:8s}: ' + ' '.join(f'k={k}: {v:.3f}' for k, v in d.items() if k in (1, 2, 3, 5, 8, 9)))
    if abl:
        base = res['all100']['ln']; b8, b16 = base[res['all100']['ks'].index(8)], base[-1]
        say(f'\nleave-one-material-out: ln RMSE of the table on the 99-material basis (all 100: k = 8: {b8:.3f}, k = 16: {b16:.3f}); largest increases:')
        for n, (r8, r16) in sorted(abl.items(), key=lambda kv: -(kv[1][1] - b16))[:10]:
            say(f'  without {n:26s} k = 8: {r8:.3f} ({r8 - b8:+.3f})   k = 16: {r16:.3f} ({r16 - b16:+.3f})')
        d16 = np.array([r[1] - b16 for r in abl.values()])
        say(f'  median change {np.median(d16):+.4f}, std {d16.std():.4f}; {int((d16 > 0.02).sum())} of {len(d16)} materials raise the k = 16 residual by more than 0.02')
    open(os.path.join(out_dir, 'report.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')


def _plot(res, ctrl, out_dir):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    style = {'all100': ('k', '-', 'o'), 'metals': ('tab:red', '-', 's'), 'metals, no aluminium': ('tab:red', '--', 's'),
             'metals+finishes': ('tab:blue', '-', '^'), 'metals+finishes, no aluminium': ('tab:blue', '--', '^')}
    for a, key, ttl in zip(axs, ('ln', 'lin'), ('ln RMSE of the log BRDF', 'relative RMSE of f cos cos')):
        for size, col in zip(ctrl, ('tab:orange', 'tab:cyan')):
            c = ctrl[size]; m, s = c[key].mean(0), c[key].std(0)
            a.fill_between(c['ks'], c[key].min(0), c[key].max(0), color=col, alpha=0.15, label=f'random {size}, min-max')
            a.plot(c['ks'], m, '-', c=col, lw=1, label=f'random {size}, mean')
        for label, r in res.items():
            col, ls, mk = style[label]
            a.plot(r['ks'], r[key], ls, c=col, marker=mk, ms=4, label=f'{label} (n = {r["n"]})')
        a.set(xlabel='k', ylabel=ttl, title=f'hpBRDF aluminium: {ttl}'); a.grid(alpha=0.3)
    axs[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, 'subsets.png'), dpi=130); plt.close(fig)
