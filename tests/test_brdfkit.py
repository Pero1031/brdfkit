"""
Self-contained tests on synthetic tables (no dataset needed). Run with  python tests/run_tests.py  (or pytest).
Tests that need a real dataset file are skipped when the file is absent.
"""
import os
import sys
import tempfile
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from brdfkit import BRDFTable                                  # noqa: E402
from brdfkit.table import WAVELENGTH, RGB                      # noqa: E402
from brdfkit import geometry as G                              # noqa: E402
from brdfkit.io import hpbrdf                                  # noqa: E402
from brdfkit.compress import compress, evaluate                # noqa: E402
from brdfkit.decode import decode, PCADecoder                  # noqa: E402
from brdfkit.export import hpbpca                              # noqa: E402
from brdfkit.analysis import lu_chipman as LC                  # noqa: E402
from brdfkit.analysis import merl_pca as MP                    # noqa: E402
from brdfkit.analysis import merl_project as MJ                # noqa: E402
from brdfkit import color                                      # noqa: E402
from brdfkit.analysis import merl_subsets as MS                # noqa: E402
from brdfkit.analysis import floor as FL                       # noqa: E402
from brdfkit.io.merl import grid_axes as merl_grid_axes        # noqa: E402

# real datasets (optional): data/kaist/hpbrdf/Aluminum.hpbrdf and the MERL .binary files in data/merl (or $BRDFKIT_MERL_DIR)
DATA_HPBRDF = os.environ.get('BRDFKIT_HPBRDF', os.path.join(ROOT, 'data', 'kaist', 'hpbrdf', 'Aluminum.hpbrdf'))
DATA_MERL_DIR = os.environ.get('BRDFKIT_MERL_DIR', os.path.join(ROOT, 'data', 'merl'))
DATA_MERL = os.path.join(DATA_MERL_DIR, 'aluminium.binary')


# ----------------------------------------------------------------------------- synthetic tables
def synthetic_polarized(n_phi=13, n_thd=10, n_thh=12, wvls=(450, 500, 550, 600, 650, 700), seed=0):
    """smooth, low-rank, physically plausible polarized table with a NaN corner"""
    phi_d = np.linspace(-np.pi, np.pi, n_phi); theta_d = np.linspace(0, np.pi / 2, n_thd); theta_h = np.linspace(0, np.pi / 2, n_thh) ** 2 / (np.pi / 2)
    wv = np.array(wvls, float)
    P, D, H = np.meshgrid(phi_d, theta_d, theta_h, indexing='ij')
    spec = 0.6 + 0.4 * (wv - wv[0]) / (wv[-1] - wv[0])                                     # per-wavelength tint
    lobe = 0.05 + 8.0 * np.exp(-(H / 0.08) ** 2) * (1 + 0.3 * np.cos(P)) * (1 + 0.5 * D)   # angular shape
    m00 = lobe[..., None] * spec                                                           # [P,D,H,C]
    Dia = 0.3 * np.sin(D) ** 2; Ret = np.pi * (1 - 0.7 * (D / (np.pi / 2)) ** 2); Dep = 0.1 * H / (np.pi / 2)
    M = np.zeros(m00.shape + (4, 4), np.float32)
    M[..., 0, 0] = 1; M[..., 0, 1] = M[..., 1, 0] = Dia[..., None]; M[..., 1, 1] = 1 - Dep[..., None]
    M[..., 2, 2] = M[..., 3, 3] = (np.cos(Ret) * (1 - Dep))[..., None]
    M[..., 2, 3] = (np.sin(Ret) * (1 - Dep))[..., None]; M[..., 3, 2] = -M[..., 2, 3]
    M *= m00[..., None, None]
    M[:, -2:, -3:] = np.nan                                                                 # unmeasured grazing corner
    return BRDFTable(phi_d, theta_d, theta_h, wv, WAVELENGTH, True, dict(source='synthetic'), values=M)


def synthetic_rgb(n_phi=11, n_thd=9, n_thh=9):
    phi_d = np.linspace(0, np.pi, n_phi); theta_d = np.linspace(0, np.pi / 2, n_thd); theta_h = np.linspace(0, np.pi / 2, n_thh)
    P, D, H = np.meshgrid(phi_d, theta_d, theta_h, indexing='ij')
    lobe = 0.02 + 3.0 * np.exp(-(H / 0.15) ** 2) * (1 + 0.2 * np.cos(P))
    vals = (lobe[..., None] * np.array([1.0, 0.8, 0.5])).astype(np.float32)
    vals[:, -1, :] = np.nan
    return BRDFTable(phi_d, theta_d, theta_h, np.arange(3.0), RGB, False, dict(source='synthetic'), values=vals)


# ----------------------------------------------------------------------------- tests
def test_tensor_file_roundtrip():
    t = synthetic_polarized()
    with tempfile.TemporaryDirectory() as d:
        fn = os.path.join(d, 'syn.hpbrdf')
        hpbrdf.write(fn, t.values, t.phi_d, t.theta_d, t.theta_h, t.channels)
        t2 = hpbrdf.load(fn)
        assert t2.n_phi == t.n_phi and t2.slice_shape == t.slice_shape
        assert np.allclose(t2.channels, t.channels) and np.allclose(t2.theta_h, t.theta_h, atol=1e-6)
        a, b = t.slice(5), t2.slice(5)
        assert np.array_equal(np.isnan(a), np.isnan(b)) and np.allclose(np.nan_to_num(a), np.nan_to_num(b))


def _compress_roundtrip(t, k1, k2, tol_m00, tol_ratio=None):
    with tempfile.TemporaryDirectory() as d:
        npz = os.path.join(d, 'c.npz')
        compress(t, npz, k1=k1, k2=k2, dtype='float32', log=lambda *a: None)
        m = evaluate(t, npz, slices=list(range(t.n_phi)), log=lambda *a: None)
        assert m['m00_rel_rmse'] < tol_m00, m
        if t.polarized:
            assert m['ratio_rmse'].max() < tol_ratio, m['ratio_rmse']
        # exact bins: decoding the stored coefficients reproduces measured bins
        with np.load(npz) as zf:
            z = dict(zf)
        i = t.n_phi // 2; sl = t.slice(i); v = z['valid'][i]
        rec = decode(z['coef_logm00'][i], z['coef_ratio'][i], z)
        o, r = sl[v], rec[v]
        assert np.allclose(t.m00(r), t.m00(o), rtol=tol_m00 * 3, atol=1e-6)
        # export -> read back -> same decode
        bin_ = os.path.join(d, 'c.hpbpca')
        hpbpca.write(npz, bin_, dtype='float32')
        zb = hpbpca.read(bin_)
        rec2 = decode(zb['coef_logm00'][i], zb['coef_ratio'][i], zb)
        assert np.allclose(rec2, rec, rtol=1e-5, atol=1e-6)
        assert zb['polarized'] == t.polarized and np.array_equal(zb['valid'], z['valid'])
        # decoder at random directions is finite and positive
        dec = PCADecoder(npz)
        rng = np.random.default_rng(0)
        def hemi(n):
            u = rng.random((n, 2)); th = np.arccos(u[:, 0]); ph = 2 * np.pi * u[:, 1]
            return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
        r = dec.eval(hemi(200), hemi(200), float(t.channels[len(t.channels) // 2]))
        assert np.isfinite(r['m00']).all() and (r['m00'] > 0).all()
        if t.polarized:
            assert np.isfinite(r['M']).all() and np.allclose(r['M'][:, 0, 0], r['m00'])
        return npz


def test_compress_polarized():
    _compress_roundtrip(synthetic_polarized(), k1=3, k2=8, tol_m00=2e-3, tol_ratio=2e-3)


def test_compress_rgb_unpolarized():
    t = synthetic_rgb()
    _compress_roundtrip(t, k1=2, k2=32, tol_m00=2e-3)


def test_geometry_rotations():
    rng = np.random.default_rng(3)
    L = G.normalize(rng.normal(size=(50, 3))); L[:, 2] = np.abs(L[:, 2]); V = G.normalize(rng.normal(size=(50, 3))); V[:, 2] = np.abs(V[:, 2])
    L = G.normalize(L); V = G.normalize(V)
    g = G.geometry(L, V)
    assert np.all(g['theta_h'] >= 0) and np.all(g['theta_h'] <= np.pi / 2 + 1e-9) and np.all(np.abs(g['phi_d']) <= np.pi + 1e-9)
    assert np.allclose((g['x_in'] * -L).sum(-1), 0, atol=1e-6) and np.allclose((g['x_out'] * V).sum(-1), 0, atol=1e-6)
    M = rng.normal(size=(50, 4, 4))
    same = G.rotate_mueller_basis(M, -L, g['x_in'], g['x_in'], V, g['x_out'], g['x_out'])
    assert np.allclose(same, M, atol=1e-9)
    sp = G.mueller_to_sp(M, g['x_in'], g['x_out'], L, V)
    assert np.allclose(sp[:, 0, 0], M[:, 0, 0]) and np.allclose(sp[:, 0, 1] ** 2 + sp[:, 0, 2] ** 2 + sp[:, 0, 3] ** 2,
                                                                 M[:, 0, 1] ** 2 + M[:, 0, 2] ** 2 + M[:, 0, 3] ** 2)


def test_lu_chipman():
    assert LC.self_test(verbose=False)
    rs, rp = LC.fresnel_rs_rp(np.cos(np.radians(60)), 1.02 + 6.63j)
    d = LC.lu_chipman(LC.mueller_from_amplitudes(rs, rp))
    assert 0 < d['D'] < 0.2 and np.pi / 2 < d['R'] < np.pi and abs(d['detmR'] - 1) < 1e-6


def test_hpbrdf_loader_real():
    if not os.path.exists(DATA_HPBRDF):
        return 'skipped (no hpBRDF file)'
    t = hpbrdf.load(DATA_HPBRDF)
    assert t.polarized and t.n_ch == 68 and t.slice_shape == (91, 91, 68, 4, 4)
    sl = t.slice(180); v = t.valid(sl)
    assert 0.5 < v.mean() < 1.0 and np.nanmax(t.m00(sl)) > 1


def synthetic_material_set(n_mat=12, n_phi=12, n_thd=10, n_thh=10, rank=3, noise=0.01, seed=0):
    """RGB tables whose log BRDF lives (up to noise) in a `rank`-dimensional affine space -> (list of (name, table), coefficients)"""
    rng = np.random.default_rng(seed)
    phi_d = np.linspace(0, np.pi, n_phi, endpoint=False); theta_d = np.linspace(0, np.pi / 2, n_thd, endpoint=False)
    theta_h = (np.arange(n_thh) / n_thh) ** 2 * (np.pi / 2)
    P, D, H = np.meshgrid(phi_d, theta_d, theta_h, indexing='ij')
    base = np.log(0.05 + 3.0 * np.exp(-(H / 0.2) ** 2))[..., None] + np.log(np.array([1.0, 0.8, 0.6]))
    B = np.stack([np.exp(-(H / 0.1) ** 2)[..., None] * np.array([1.0, 1.0, 1.0]),
                  (D / (np.pi / 2))[..., None] * np.array([1.0, 0.5, 0.2]),
                  (np.cos(P) * 0.3)[..., None] * np.array([0.2, 1.0, 0.5])])[:rank]
    A = rng.normal(size=(n_mat, rank)) * np.array([2.0, 1.0, 0.5])[:rank]
    mats = []
    for i in range(n_mat):
        x = base + np.tensordot(A[i], B, 1) + noise * rng.normal(size=base.shape)
        v = np.exp(x).astype(np.float32); v[:, -1, -2:] = np.nan
        mats.append((f'm{i:02d}', BRDFTable(phi_d, theta_d, theta_h, np.arange(3.0), RGB, False, dict(source='synthetic'), values=v)))
    return mats, A


def test_merl_pca_cos_weights():
    """cos_i cos_o from the half/diff inverse agrees with brdfkit.geometry on the forward conversion"""
    rng = np.random.default_rng(5)
    pd, td, th = rng.uniform(0, np.pi, 200), rng.uniform(0, np.pi / 2, 200), rng.uniform(0, np.pi / 2, 200)
    w = np.array([MP.cos_weights([p], [d], [h])[0, 0, 0] for p, d, h in zip(pd, td, th)])
    d = np.stack([np.sin(td) * np.cos(pd), np.sin(td) * np.sin(pd), np.cos(td)], 1)
    wi = np.stack([d[:, 0] * np.cos(th) + d[:, 2] * np.sin(th), d[:, 1], -d[:, 0] * np.sin(th) + d[:, 2] * np.cos(th)], 1)
    hv = np.stack([np.sin(th), 0 * th, np.cos(th)], 1); wo = 2 * (wi * hv).sum(1, keepdims=True) * hv - wi
    ok = (wi[:, 2] > 1e-3) & (wo[:, 2] > 1e-3)
    assert ok.sum() > 50 and np.allclose(w[ok], wi[ok, 2] * wo[ok, 2]) and np.all(w[~ok] == 0)
    g = G.geometry(wi[ok], wo[ok])
    assert np.allclose(g['theta_h'], th[ok], atol=1e-6) and np.allclose(g['theta_d'], td[ok], atol=1e-6)


def _check_material_pca(mapping):
    mats, A = synthetic_material_set()
    model, Xc = MP.fit(mats, stride=1, mapping=mapping, k=64, log=lambda *a: None)
    n = len(mats)
    assert model.coef.shape == (n, n) and model.basis.shape[1] == n - 2 and model.loo_rmse.shape == (n, n - 1)
    assert model.evr[:3].sum() > 0.97, model.evr[:5]
    assert np.median(model.recon_rmse[:, 3]) < 0.02 and np.median(model.loo_rmse[:, 3]) < 0.05, (model.recon_rmse[:, 3], model.loo_rmse[:, 3])
    assert np.all(np.diff(model.recon_rmse, axis=1) <= 1e-9) and np.all(model.recon_rmse[:, n - 1] < 1e-4)
    # the 3-d latent space is recovered: coefficients are an affine function of the generating A
    Z = np.c_[A, np.ones(n)]; fitc = np.linalg.lstsq(Z, model.coef[:, :3], rcond=None)[0]
    assert np.abs(Z @ fitc - model.coef[:, :3]).max() < 0.2 * np.abs(model.coef[:, :3]).max()
    # project() of a fitted table reproduces its stored coefficients; reconstruct() reproduces the table
    i = 4; v = mats[i][1].values
    assert np.allclose(model.project(v, k=5), model.coef[i, :5], atol=1e-3)
    rec = model.reconstruct(model.coef[i], k=n - 2)
    m = model.valid
    assert np.array_equal(np.isnan(rec).any(-1), ~m)
    assert np.allclose(rec[m], v[m], rtol=0.05 if mapping == 'nielsen' else 0.02, atol=1e-3)
    rec3 = model.reconstruct(model.coef[i], k=3)
    assert np.sqrt(np.mean((np.log(rec3[m]) - np.log(v[m])) ** 2)) < 0.1
    # save / load round trip
    with tempfile.TemporaryDirectory() as d:
        fn = model.save(os.path.join(d, 'm.npz')); m2 = MP.MaterialPCA.load(fn)
        assert m2.names == model.names and m2.mapping == mapping and m2.stride == (1, 1, 1)
        assert np.allclose(m2.project(v, k=5), model.coef[i, :5], atol=1e-3)
        assert np.allclose(m2.reconstruct(model.coef[i], k=3)[m], rec3[m], rtol=1e-4)
    return model


def test_merl_pca_synthetic_log():
    _check_material_pca('log')


def test_merl_pca_synthetic_nielsen():
    m = _check_material_pca('nielsen')
    assert m.ref is not None and m.ref.shape == (m.dim,)


def test_merl_pca_real():
    if not os.path.isdir(DATA_MERL_DIR):
        return 'skipped (no MERL directory)'
    files = [os.path.join(DATA_MERL_DIR, f) for f in ('aluminium.binary', 'alum-bronze.binary', 'black-obsidian.binary',
                                                      'beige-fabric.binary', 'gold-metallic-paint.binary')]
    files = [f for f in files if os.path.exists(f)]
    if len(files) < 4:
        return 'skipped (too few MERL files)'
    with tempfile.TemporaryDirectory() as d:
        model = MP.run_report(files, d, stride=6, k=2, plot=False, log=lambda *a: None)
        assert os.path.exists(os.path.join(d, 'merl_pca.npz')) and os.path.exists(os.path.join(d, 'merl_pca_report.txt'))
        assert model.grid_shape == (30, 15, 15) and model.dim == model.valid.sum() * 3 and 0.6 < model.valid.mean() < 0.9
        assert model.evr[0] > 0.3 and np.all(model.recon_rmse[:, len(files) - 1] < 1e-3)
        from brdfkit.io import merl
        i = model.names.index('aluminium')                                       # run_report sorts materials by name
        assert np.allclose(model.project(merl.load(files[0]).values), model.coef[i, :2], atol=1e-2)
    return f'{len(files)} materials, evr1 {model.evr[0]:.2f}'


def test_merl_loader_real():
    if not os.path.exists(DATA_MERL):
        return 'skipped (no MERL file)'
    from brdfkit.io import merl
    t = merl.load(DATA_MERL)
    assert not t.polarized and t.channel_kind == RGB and t.values.shape == (180, 90, 90, 3)
    sl = t.slice(90); v = t.valid(sl)
    assert v.mean() > 0.5 and np.nanmax(sl) > 1 and np.nanmin(sl[v]) >= 0


# ----------------------------------------------------------------------------- hpBRDF -> MERL PCA space
def test_color_weights():
    """band weights: columns sum to 1, flat spectrum -> equal RGB, constants identical to the renderer's MeasuredHpBRDF.hpp"""
    import re
    wv = np.arange(414, 951, 8.0)
    w = color.rgb_band_weights(wv)
    assert w.shape == (68, 3) and np.allclose(w.sum(0), 1.0)
    assert np.allclose(np.ones(68) @ w, 1.0)
    red = np.where(wv > 600, 1.0, 0.0) @ w
    blue = np.where(wv < 480, 1.0, 0.0) @ w
    assert red[0] > 0.3 and red[0] > red[2] and blue[2] > 0.3 and blue[2] > blue[0], (red, blue)
    hdr = open(os.path.join(ROOT, 'integrations', 'gold_rayt', 'include', 'Materials', 'MeasuredHpBRDF.hpp'), encoding='utf-8').read()
    for axis in 'xyz':
        line = re.search(r'(?m)^\s*' + axis + r' = ([^;]*);', hdr[hdr.index('static void cmf'):]).group(1)
        lobes = [(-1.0 if sgn.strip() == '-' else 1.0, float(a), float(mu), float(s1), float(s2))
                 for sgn, a, mu, s1, s2 in re.findall(r'([+-]?)\s*(\d+\.\d+) \* lobe\(lam, ([\d.]+), ([\d.]+), ([\d.]+)\)', line)]
        assert len(lobes) == len(color.CMF_LOBES[axis]), (axis, line)
        for (sg, a, mu, s1, s2), (a0, mu0, s10, s20) in zip(lobes, color.CMF_LOBES[axis]):
            assert abs(sg * a - a0) < 1e-12 and (mu, s1, s2) == (mu0, s10, s20), (axis, line)
    for v in color.XYZ_TO_SRGB.ravel():
        assert f'{abs(v):.4f}' in hdr, v


def _merl_like_table(asym=0.0, nan_corner=True):
    """unpolarized 4-band table on a grid containing the MERL grid (phi_d every 1 deg over [-180, 180], theta_d and the squared
    theta_h axis up to and including 90 deg); f obeys the symmetries f(phi) = f(-phi) = f(phi + pi) up to the antisymmetric
    term asym * sin(phi), which the phi-equivalent average removes"""
    from brdfkit.io import merl
    pm, tdm, thm = merl.grid_axes()
    phi = np.radians(np.arange(-180, 181.0)); td = np.append(tdm, np.pi / 2); th = np.append(thm, np.pi / 2)
    wv = np.array([450.0, 550.0, 600.0, 650.0])

    def f(j):
        D, H = np.meshgrid(td, th, indexing='ij')
        v = (0.05 + 4.0 * np.exp(-(H / 0.1) ** 2) * (1 + 0.5 * D)) * (1 + 0.3 * np.cos(2 * phi[j])) + asym * np.sin(phi[j])
        v = np.repeat(v[..., None], 4, -1).astype(np.float32)
        if nan_corner:
            v[-3:, -3:] = np.nan
        return v
    return BRDFTable(phi, td, th, wv, WAVELENGTH, False, dict(source='synthetic'), slice_fn=f), (pm, tdm, thm)


def test_merl_grid_rgb():
    t, (pm, tdm, thm) = _merl_like_table(asym=0.2)
    full, info = MJ.merl_grid_rgb(t, 'full', log=lambda *a: None)
    none, _ = MJ.merl_grid_rgb(t, 'none', log=lambda *a: None)
    assert full.values.shape == (180, 90, 90, 3) and none.values.shape == full.values.shape
    D, H = np.meshgrid(tdm, thm, indexing='ij')
    for b in (0, 17, 90, 133, 179):
        g = (0.05 + 4.0 * np.exp(-(H / 0.1) ** 2) * (1 + 0.5 * D)) * (1 + 0.3 * np.cos(2 * pm[b]))
        ok = np.isfinite(full.values[b]).all(-1)
        assert ok.sum() > 0.9 * ok.size
        # flat spectrum -> equal RGB; the symmetric average removes the antisymmetric term, 'none' keeps it
        assert np.allclose(full.values[b][ok], np.repeat(g[ok][:, None], 3, 1), rtol=2e-5, atol=2e-5)
        assert np.allclose(none.values[b][ok], np.maximum(np.repeat(g[ok][:, None], 3, 1) + 0.2 * np.sin(pm[b]), 0), rtol=2e-5, atol=2e-5)
    assert np.isnan(full.values[:, -2:, -2:]).all() and np.isfinite(full.values[:, :-2]).all()    # the NaN corner, cut to the MERL range
    assert np.array_equal(full.values[10], full.values[170], equal_nan=True) and np.isfinite(info['spread']).any() and np.nanmax(info['spread']) > 0.01
    # a table whose axes do not contain the MERL grid is rejected
    bad = BRDFTable(t.phi_d, t.theta_d, t.theta_h * 0.99, t.channels, WAVELENGTH, False, {}, slice_fn=t.slice_fn)
    try:
        MJ.merl_grid_rgb(bad, log=lambda *a: None)
    except ValueError as e:
        assert 'theta_h' in str(e)
    else:
        raise AssertionError('mismatching axis accepted')


def test_project_masked():
    mats, A = synthetic_material_set()
    model, Xc = MP.fit(mats, stride=1, mapping='log', k=6, log=lambda *a: None)
    v = mats[3][1].values
    c, info = model.project_masked(v, 4)
    assert np.allclose(c, model.project(v, 4), atol=1e-3) and info['n_used'] == info['n_total'] and info['cond'] < 1.001
    # unmeasured bins inside the model mask: the masked fit stays close, the plain fit refuses
    w = v.copy(); w[:, 2, 2:5] = np.nan
    c2, info2 = model.project_masked(w, 4)
    assert info2['n_used'] < info2['n_total'] and np.abs(c2 - c).max() < 0.1 * np.abs(c).max(), (c, c2)
    try:
        model.project(w)
    except ValueError:
        pass
    else:
        raise AssertionError('project accepted unmeasured bins')


def test_merl_project_real():
    """the hpBRDF axes contain the MERL grid, and one slice of the real table gives a sane RGB slice"""
    if not os.path.exists(DATA_HPBRDF):
        return 'skipped (no hpBRDF file)'
    from brdfkit.io import merl
    t = hpbrdf.load(DATA_HPBRDF)
    pm, tdm, thm = merl.grid_axes()
    MJ._index_of(t.theta_d, tdm, 'theta_d'); MJ._index_of(t.theta_h, thm, 'theta_h')
    j = int(np.argmin(np.abs(t.phi_d - np.pi / 2)))
    rgb = t.m00(t.slice(j)) @ color.rgb_band_weights(t.channels)
    fin = np.isfinite(rgb).all(-1)
    assert fin.mean() > 0.9 and rgb[fin].max() > 1 and np.median(rgb[fin]) < 1
    pos = rgb[fin][(rgb[fin] > 1e-3).all(-1)]
    peak = pos[pos.mean(-1) > np.percentile(pos.mean(-1), 95)]                          # the specular region of aluminium is nearly grey
    assert np.abs(np.median(np.log(peak[:, 0] / peak[:, 2]))) < 0.3


def test_subset_residuals():
    """residuals() of a table that lives in the basis: tiny at the true rank, decreasing in k, measured on the common mask only"""
    mats, A = synthetic_material_set(n_mat=14)
    model, _ = MP.fit(mats[:10], stride=1, mapping='log', k=8, log=lambda *a: None)
    hp = mats[12][1].values.copy()                                                  # a material the basis has not seen
    common = model.valid & np.isfinite(hp).all(-1)
    ln, lin, zmax, zc = MS.residuals(model, hp, common, [1, 2, 3, 5])
    assert np.all(np.diff(ln) <= 1e-6) and ln[2] < 0.1 and lin[2] < 0.1 and ln[0] > 3 * ln[2], (ln, lin)
    assert np.isfinite(zmax) and 1 <= zc <= 5
    # unmeasured bins in the table are skipped, not propagated
    hp[3, 2:5, 2:5] = np.nan
    ln2 = MS.residuals(model, hp, model.valid & np.isfinite(hp).all(-1), [3])[0]
    assert np.isfinite(ln2).all() and abs(ln2[0] - ln[2]) < 0.05
    for name in MS.METALS + MS.FINISHES:
        assert isinstance(name, str)


def test_grid_stats():
    """peak shape and floor statistics of a table on the MERL grid"""
    th = np.arange(90)
    core = np.where(th < 6, 10.0, 10.0 * np.exp(-(th - 5) / 3.0))                      # flat for the first 6 bins, then decays
    v = np.full((180, 90, 90, 3), 0.01, np.float32)
    v[:, :, :, 1] = np.maximum(core, 0.01)[None, None, :]
    v[:, :, 60:, 1] = 0.002
    v[10, 5, 3, 1] = 50.0
    st = MJ.grid_stats(v)
    assert st['peak_max'] == 50.0 and st['peak_core'] == 10.0 and abs(st['flat_ratio'] - 1.0) < 1e-6 and abs(st['floor'] - 0.002) < 1e-9
    v[:, :, :, 1] = np.maximum(10.0 * np.exp(-th / 3.0), 0.01)[None, None, :]
    assert MJ.grid_stats(v)['flat_ratio'] < 0.2


def _floored_table():
    """black-glass-like polarized table: a narrow lobe that is zero far away, plus an unpolarized offset c(lambda) on m00"""
    phi_d = np.linspace(-np.pi, np.pi, 13); theta_d = np.radians(np.arange(0, 90, 5.0)); theta_h = (np.arange(30) / 30) ** 2 * (np.pi / 2)
    wv = np.array([450.0, 550.0, 650.0]); c = np.array([0.02, 0.006, 0.004])
    P, D, Hh = np.meshgrid(phi_d, theta_d, theta_h, indexing='ij')
    true = 5.0 * np.exp(-(Hh / 0.05) ** 2)[..., None] * np.ones(3)                          # [P, D, H, C]
    M = np.zeros(true.shape + (4, 4), np.float32)
    for k, d in enumerate((1.0, 0.9, -0.9, -0.9)):
        M[..., k, k] = true * d
    M[..., 0, 0] += c
    return BRDFTable(phi_d, theta_d, theta_h, wv, WAVELENGTH, True, dict(source='synthetic'), values=M), true, c


def test_floor_estimate_subtract():
    """the floor is estimated far from the lobe, is unpolarized, and subtracting it restores m00 (clamped where it would vanish)"""
    t, true, c = _floored_table()
    fl = FL.estimate(t, phi_step=1, log=lambda *a: None)
    assert np.allclose(fl['offset'], c, rtol=1e-4) and fl['norm_absmed'][:, 1:, 1:].max() < 1e-6
    s = FL.subtract(t, fl)
    m = s.values[..., 0, 0]
    big = true > 10 * c
    assert np.allclose(m[big], true[big], rtol=1e-4)                                          # lobe: offset removed
    far = true < 1e-12
    assert np.allclose(m[far], (FL.KEEP * c * np.ones_like(true))[far], rtol=1e-4)            # far: clamped at keep * measured
    assert np.allclose(s.values[big][..., 1:, :], t.values[big][..., 1:, :])                 # in the lobe only m00 changes
    row = np.sqrt((s.values[..., 0, 1:] ** 2).sum(-1)); assert (m >= row - 1e-7).all()        # diattenuation stays <= 1
    raw = FL.subtract(t, fl, max_gain=None)
    assert np.allclose(raw.values[..., 0, 0], m) and np.array_equal(raw.values[..., 1:, :], t.values[..., 1:, :])
    x = t.values.copy(); x[..., 1, 1] = 0.5 * x[..., 0, 0]                                     # weak signal under the offset
    tx = BRDFTable(t.phi_d, t.theta_d, t.theta_h, t.channels, WAVELENGTH, True, {}, values=x)
    g = FL.subtract(tx, fl, max_gain=4.0).values
    with np.errstate(all='ignore'):
        rise = (g[..., 1, 1] / g[..., 0, 0]) / (x[..., 1, 1] / x[..., 0, 0])
    assert np.nanmax(rise) <= 4.0 + 1e-4 and np.nanmax(np.abs(g[..., 1, 1] / g[..., 0, 0])) <= 1 + 1e-6
    assert np.allclose(g[big][..., 0, 0], true[big], rtol=1e-4)                               # the lobe is untouched by the cap
    lazy = FL.subtract(BRDFTable(t.phi_d, t.theta_d, t.theta_h, t.channels, WAVELENGTH, True, {}, slice_fn=lambda i: t.values[i]), fl)
    assert np.allclose(lazy.slice(3), s.values[3])
    bad = dict(fl, channels=np.array([400.0, 500.0, 600.0]))
    try:
        FL.subtract(t, bad); raise AssertionError('channel mismatch not detected')
    except ValueError:
        pass
    with tempfile.TemporaryDirectory() as d:
        FL.save(os.path.join(d, 'f.npz'), fl)
        assert np.allclose(FL.subtract(t, os.path.join(d, 'f.npz')).values, s.values)


def test_floor_rgb():
    """floor removal on a MERL-grid RGB table: subtract the offset taken to RGB, clamp at keep * value, keep NaN"""
    v = np.full((180, 90, 90, 3), 0.01, np.float32); v[:, :, :5] = 3.0; v[0, 0, 0] = np.nan
    t = BRDFTable(*merl_grid_axes(), np.arange(3.0), RGB, False, {}, values=v)
    fl = dict(offset=np.array([0.008, 0.009, 0.0095]), channels=np.arange(3.0), channel_kind=RGB)
    r = MJ.floor_rgb(t, fl).values
    assert np.isnan(r[0, 0, 0]).all() and np.allclose(r[1, 1, 1], 3.0 - fl['offset'])
    assert np.allclose(r[1, 1, 50], 0.01 - fl['offset'], atol=1e-7)
    v[:, :, 50:] = 0.005
    r = MJ.floor_rgb(t, fl, keep=0.02).values
    assert np.allclose(r[1, 1, 60], 0.0001)
