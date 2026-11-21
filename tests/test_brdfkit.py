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
from brdfkit.analysis import factor as FA                      # noqa: E402
from brdfkit.analysis import latent_edit as LE                 # noqa: E402
from brdfkit.analysis import calibrate as CA                  # noqa: E402
from brdfkit.analysis import swap as SW                       # noqa: E402
from brdfkit.analysis import gloss as GL                      # noqa: E402
from brdfkit.analysis import quirks as QK                     # noqa: E402
from brdfkit.analysis import slots as SL                      # noqa: E402
from brdfkit.analysis import clean as CL                      # noqa: E402
from brdfkit.analysis import umov as UM                       # noqa: E402
from brdfkit.analysis import predict as PR                    # noqa: E402
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


def _separable_table(noise=0.0, seed=0):
    """m00 = A(angles) S(theta_d, lambda), m_ij/m00 = N(theta_d, lambda): exactly the factor model (plus optional noise)"""
    rng = np.random.default_rng(seed)
    phi_d = np.linspace(-np.pi, np.pi, 9); theta_d = np.radians(np.arange(0, 91, 5.0)); theta_h = (np.arange(16) / 15) ** 2 * (np.pi / 2)
    wv = np.linspace(420, 900, 7)
    P, D, Hh = np.meshgrid(phi_d, theta_d, theta_h, indexing='ij')
    A = 0.02 + 6.0 * np.exp(-(Hh / 0.1) ** 2) * (1 + 0.2 * np.cos(P))
    td = theta_d[:, None]
    S = 0.5 + 0.3 * np.sin(td) * (wv[None] - 420) / 480 + 0.1 * np.cos(wv[None] / 100)
    N = np.zeros((len(theta_d), len(wv), 15))
    N[..., 0] = 0.4 * np.sin(td) ** 2 * (wv[None] / 900); N[..., 4] = 0.95; N[..., 9] = -0.9 * np.cos(2 * td); N[..., 14] = N[..., 9]
    m00 = A[..., None] * S[None, :, None, :]
    M = np.empty(m00.shape + (16,), np.float32)
    M[..., 0] = m00 * (1 + noise * rng.standard_normal(m00.shape)); M[..., 1:] = m00[..., None] * N[None, :, None]
    M = M.reshape(m00.shape + (4, 4)); M[0, -1, -1] = np.nan
    return BRDFTable(phi_d, theta_d, theta_h, wv, WAVELENGTH, True, dict(source='synthetic'), values=M), A, S, N


def test_factor_exact():
    """a separable table is recovered: A x S reproduces m00, N matches, the view rebuilds the table, save/load/io.load round-trip"""
    t, A, S, N = _separable_table()
    m0, m5 = FA.fit(t, measured_theta_d=(0, 90), smooth_deg=(0.0, 10.0), log=lambda *a: None)
    assert m0['errors']['m00_rel_rmse']['all_theta_d']['all'] < 1e-5 and m0['errors']['N_rmse']['all_theta_d']['all'] < 1e-5
    assert np.allclose(m0['N'], N, atol=1e-5)
    yw = FA.luminance_weights(t.channels); assert np.allclose(m0['S'] @ yw, 1.0)
    assert np.allclose(m0['S'], S / (S @ yw)[:, None], rtol=1e-5)                         # S up to the luminance normalization
    assert m5['errors']['m00_rel_rmse']['all_theta_d']['all'] > m0['errors']['m00_rel_rmse']['all_theta_d']['all']
    v = FA.view(m0)
    ok = np.isfinite(t.values)
    assert np.allclose(v.slice(4)[ok[4]], t.values[4][ok[4]], rtol=1e-4, atol=1e-6) and np.isnan(v.slice(0)[-1, -1]).all()
    from brdfkit.io import load as io_load
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'm.npz'); FA.save(p, m0)
        assert FA.is_model(p) and np.allclose(io_load(p).slice(2), v.slice(2), rtol=2e-3, atol=1e-5)   # A stored as float16


def test_factor_measured_range():
    """outside the measured theta_d range S and N are copied from the nearest measured theta_d"""
    t, A, S, N = _separable_table(noise=0.01)
    m = FA.fit(t, measured_theta_d=(20, 80), log=lambda *a: None)[0]
    td = np.degrees(t.theta_d); lo = np.nonzero(td >= 20)[0][0]; hi = np.nonzero(td <= 80)[0][-1]
    assert np.allclose(m['S'][:lo], m['S'][lo]) and np.allclose(m['S'][hi + 1:], m['S'][hi]) and np.allclose(m['N'][:lo], m['N'][lo])
    e = m['errors']['m00_rel_rmse']
    assert e['measured_theta_d']['all'] < 0.012 and e['all_theta_d']['all'] > 1.5 * e['measured_theta_d']['all']   # copied rows miss the true S


def _edit_setup():
    """small log PCA of synthetic RGB materials, and a factor model whose grid contains the PCA grid (phi_d over [-180, 180])"""
    mats, _ = synthetic_material_set()
    pca, _ = MP.fit(mats, stride=1, mapping='log', k=8, log=lambda *a: None)
    phi = np.radians(np.arange(-180, 181, 15.0)); td = pca.theta_d; th = pca.theta_h
    P, D, Hh = np.meshgrid(phi, td, th, indexing='ij')
    A = (0.05 + 3.0 * np.exp(-(Hh / 0.2) ** 2)).astype(np.float32)
    S = np.ones((len(td), 3)); N = np.zeros((len(td), 3, 15))
    m = dict(kind=FA.KIND, A=A, S=S, N=N, phi_d=phi, theta_d=td, theta_h=th, channels=np.array([450.0, 550.0, 650.0]),
             measured_theta_d=np.array([0.0, 90.0]), smooth_deg=0.0, errors={}, meta={})
    return pca, m


def test_latent_edit():
    """no shift = no change; a shift along one component multiplies A by exp(that component) at the PCA grid points;
    keep_albedo restores the albedo; the grey A projects close to the base material of the synthetic set"""
    pca, m = _edit_setup()
    same, info = LE.edit(m, pca)
    assert np.allclose(same['A'], m['A']) and len(info['z']) == 8
    e, info = LE.edit(m, pca, pcs=[(2, 1.0)])
    sig = pca.coef[:, :8].std(0)
    expect = (pca.basis[:, 1] * sig[1]).reshape(-1, 3).mean(1)
    j15 = int(np.argmin(np.abs(m['phi_d'] - pca.phi_d[1])))                    # phi_d grid point of the PCA inside A's grid
    got = np.log(e['A'][j15] / m['A'][j15])
    v = pca.valid[1]
    full = np.zeros(pca.grid_shape); full[pca.valid] = expect
    assert np.allclose(got[v], full[1][v], atol=1e-4)
    k, _ = LE.edit(m, pca, pcs=[(1, 1.5)], keep_albedo_deg=30)
    assert abs(LE.albedo(k, 30) - LE.albedo(m, 30)) < 1e-6 * max(1, LE.albedo(m, 30))
    t0, _ = LE.edit(m, pca, toward=[(pca.names[3], 0.0)])
    assert np.allclose(t0['A'], m['A'])
    try:
        LE.edit(m, pca, pcs=[(9, 1.0)]); raise AssertionError('component outside k accepted')
    except ValueError:
        pass


def test_symmetrize_phi():
    """the four phi_d-equivalent bins become equal; 'median' removes a one-sided arm, 'mean' spreads it over the four;
    a symmetric lobe is unchanged; edit() symmetrizes by default, so an edit with no shift is the symmetrized model"""
    phi = np.radians(np.arange(-180, 181, 1.0))                                  # the hpBRDF phi_d axis
    idx = FA.phi_equivalents(phi)
    j = int(np.argmin(abs(phi - np.radians(-30))))
    assert sorted(np.round(np.degrees(phi[idx[j]])).astype(int)) == [-150, -30, 30, 150]
    assert set(idx[0]) == set(idx[-1])                                           # -180 and +180 deg are one angle
    base = 1.0 + 0.1 * np.cos(2 * phi)                                           # invariant under all four maps
    X = np.broadcast_to(base[:, None, None], (len(phi), 3, 4)).astype(np.float32)
    for mode in ('median', 'mean'):
        assert np.allclose(FA.symmetrize_phi(X, phi, mode), X, atol=1e-6)
    Y = X.copy(); Y[j] *= 4.0                                                    # one-sided arm at -30 deg
    med, mean = FA.symmetrize_phi(Y, phi, 'median'), FA.symmetrize_phi(Y, phi, 'mean')
    for Z in (med, mean):
        assert all(np.allclose(Z[idx[:, 0]], Z[idx[:, k]]) for k in range(1, 4))
    assert np.allclose(med, X, atol=1e-6)                                        # the arm is gone
    assert np.allclose(mean[idx[j]], X[j] * 1.75, rtol=1e-5)                     # (4 + 1 + 1 + 1) / 4 at all four angles
    try:
        FA.symmetrize_phi(X[:181], phi[180:], 'median'); raise AssertionError('half-range phi_d accepted')
    except ValueError:
        pass
    pca, m = _edit_setup()
    A = m['A'].copy(); A[int(np.argmin(abs(m['phi_d'] - np.radians(-30))))] *= 3.0
    m2 = dict(m, A=A)
    e, info = LE.edit(m2, pca)
    assert np.allclose(e['A'], m['A'], rtol=1e-5) and info['phi_symmetry'] == 'median'
    e0, _ = LE.edit(m2, pca, sym='none')
    assert np.allclose(e0['A'], A)


def _fresnel_model(nk_or_n, K, width_deg, wv=(450.0, 550.0, 650.0, 750.0)):
    """factor model of a specular reflector (dielectric n per channel, or a metal n, k .yml) on the hpBRDF grid, whose lobe
    energy is K(theta_d) times its Fresnel reflectance and whose N is the Fresnel Mueller matrix in the table frames"""
    td = np.radians(np.arange(0, 91, 1.0)); th = (np.arange(91) / 90.0) ** 2 * (np.pi / 2); phi = np.radians(np.arange(-180, 181, 15.0))
    lam = np.array(wv)
    if isinstance(nk_or_n, str):
        F = CA.metal_reflectance(td, lam, nk_or_n); met = LC.Metal(nk_or_n)
        Msp = np.stack([LC.mueller_from_amplitudes(*LC.fresnel_rs_rp(np.cos(td), met.eta(l))) for l in lam], 1)
    else:
        F, Msp = CA.dielectric(td, np.asarray(nk_or_n, float))
    yw = FA.luminance_weights(lam); FY = F @ yw
    D = np.exp(-0.5 * (np.degrees(th) / width_deg) ** 2)
    sel = np.degrees(th) <= 10.0
    D /= CA.trapezoid(D[sel] * np.cos(th[sel]) * np.sin(th[sel]), th[sel]) * 2 * np.pi             # the estimator's normalization
    with np.errstate(divide='ignore', invalid='ignore'):
        Atd = FY * K(td) / (4 * np.cos(td) ** 2)
    Atd[-1] = Atd[-2]
    A = np.broadcast_to((Atd[:, None] * D[None, :])[None], (len(phi), len(td), len(th))).astype(np.float32)
    N = CA.sp_to_table(Msp, td).reshape(len(td), len(lam), 16)[..., 1:]
    N[0] = N[1]; N[-1] = N[-2]
    return dict(kind=FA.KIND, A=A, S=F / FY[:, None], N=N, phi_d=phi, theta_d=td, theta_h=th, channels=lam,
                measured_theta_d=np.array([20.0, 80.0]), smooth_deg=0.0, errors={}, meta={'source': 'synthetic'})


def test_calibrate():
    """table <-> (s, p) frames round-trip; a synthetic glass with a known index and radiometric error K(theta_d) gives back n
    and K; applying K to a metal with the same error brings its albedo to its Fresnel reflectance"""
    td = np.radians(np.arange(0, 91, 1.0))
    F, Msp = CA.dielectric(td, np.array([1.5, 1.6]))
    back = CA.table_to_sp(CA.sp_to_table(Msp, td), td)
    ok = np.isfinite(back).all((-1, -2, -3))
    assert ok.sum() == 89 and np.allclose(back[ok], Msp[ok], atol=1e-10)
    K = lambda t: 0.5 * np.cos(t) ** 2 + 0.05
    n_true = 1.50 + 0.02 * (550.0 / np.array([450.0, 550.0, 650.0, 750.0])) ** 2
    glass = _fresnel_model(n_true, K, width_deg=0.8)
    cal = CA.estimate(glass, cutoff_deg=10, smooth_deg=0)
    assert np.allclose(cal['n'], n_true, atol=1e-3) and cal['n_rms'].max() < 1e-3
    meas = (np.degrees(td) >= 20) & (np.degrees(td) <= 80)
    assert np.allclose(cal['K'][meas], K(td[meas]), rtol=1e-4)                  # n on a 0.0002 grid, A in float32
    assert np.allclose(cal['spectral_error'][meas], 1.0, atol=1e-3) and np.nanmax(np.abs(cal['D_meas'] - cal['D_fresnel'])) < 1e-3
    with tempfile.TemporaryDirectory() as d:
        nk = os.path.join(d, 'metal.yml')
        open(nk, 'w').write('DATA:\n  - type: tabulated nk\n    data: |\n        0.40 0.50 4.8\n        0.80 2.00 8.0\n')
        metal = _fresnel_model(nk, K, width_deg=3.0)
        rows = CA.check(metal, cal, nk, angles=(30, 50), n_mc=300_000)
        for a, b, Fm, r0, r1 in rows:
            assert r0 < 0.6 and abs(r1 - 1) < 0.08, (a, r0, r1)
        p = os.path.join(d, 'cal.npz'); CA.save(p, cal)
        kept = CA.apply(metal, CA.load(p), keep_level_deg=30)
        j = int(np.argmin(abs(np.degrees(td) - 30)))
        assert np.allclose(kept['A'][:, j], metal['A'][:, j], rtol=1e-5) and kept['meta']['radiometric_calibration']['keep_level_deg'] == 30


def test_swap():
    """a full swap keeps the shape of X and takes the reflectance level, S and N of Y; a colour-only swap leaves A and N;
    mixed calibration states are refused; a Fresnel source on a calibrated shape gives that reflectance"""
    K = lambda t: 0.5 * np.cos(t) ** 2 + 0.05
    with tempfile.TemporaryDirectory() as d:
        nk = os.path.join(d, 'metal.yml')
        open(nk, 'w').write('DATA:\n  - type: tabulated nk\n    data: |\n        0.40 0.50 4.8\n        0.80 2.00 8.0\n')
        metal = _fresnel_model(nk, K, width_deg=3.0)
        glass = _fresnel_model(np.full(4, 1.5), K, width_deg=0.8)
        td = metal['theta_d']; meas = (np.degrees(td) >= 20) & (np.degrees(td) <= 80)
        out, info = SW.swap(metal, glass, cutoff_deg=10, smooth_deg=0)
        assert np.allclose(CA.lobe_energy(out, 10)[meas], CA.lobe_energy(glass, 10)[meas], rtol=1e-4)
        near = np.degrees(metal['theta_h']) < 10
        r = out['A'][:, meas][..., near] / metal['A'][:, meas][..., near]
        assert np.allclose(r, r[:1, :, :1], rtol=1e-5)                                   # same shape, rescaled per theta_d
        assert np.allclose(out['S'], glass['S']) and np.allclose(out['N'], glass['N']) and info['parts'] == ['S', 'N', 'level']
        col, _ = SW.swap(metal, glass, parts=['S'])
        assert np.allclose(col['A'], metal['A']) and np.allclose(col['N'], metal['N']) and np.allclose(col['S'], glass['S'])
        metal_c = CA.apply(metal, CA.estimate(glass, 10, 0))
        try:
            SW.swap(metal_c, glass); raise AssertionError('mixed calibration states accepted')
        except ValueError:
            pass
        fr = SW.fresnel_source(metal_c, nk)
        assert np.allclose(fr['N'][meas], metal['N'][meas], atol=1e-6) and np.allclose(fr['S'], metal['S'])
        phys, _ = SW.swap(metal_c, fr, cutoff_deg=10, smooth_deg=0)
        F = CA.metal_reflectance(td, metal['channels'], nk) @ FA.luminance_weights(metal['channels'])
        assert np.allclose(CA.lobe_energy(phys, 10)[meas], F[meas], rtol=1e-3)


def _microfacet_model(alpha, kind='beckmann'):
    """factor model whose rows are a normalized GGX or Beckmann D(theta_h) / (4 cos_i cos_o) on the hpBRDF grid (phi_d coarse)"""
    phi = np.radians(np.arange(-180, 181, 15.0)); td = np.radians(np.arange(0, 91, 1.0)); th = (np.arange(91) / 90.0) ** 2 * (np.pi / 2)
    ci, co = GL.cosines(phi, td, th)
    t2 = np.tan(th) ** 2; c4 = np.cos(th) ** 4
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        D = np.exp(-t2 / alpha ** 2) / (np.pi * alpha ** 2 * c4) if kind == 'beckmann' else alpha ** 2 / (np.pi * c4 * (alpha ** 2 + t2) ** 2)
        A = np.where((ci > 0) & (co > 0), np.broadcast_to(D, ci.shape) / (4 * ci * co), np.nan)
    lam = np.array([450.0, 550.0, 650.0])
    return dict(kind=FA.KIND, A=A.astype(np.float32), S=np.ones((len(td), 3)), N=np.zeros((len(td), 3, 15)), phi_d=phi, theta_d=td,
                theta_h=th, channels=lam, measured_theta_d=np.array([20.0, 80.0]), smooth_deg=0.0, errors={}, meta={'source': 'synthetic'})


def test_gloss():
    """Beckmann (*) Beckmann is Beckmann with the alphas added in quadrature; the slope mass is kept; a target width is hit;
    a narrower target is refused; S and N are untouched; the GGX hwhm <-> alpha conversion round-trips"""
    for a in (0.02, 0.1, 0.4):
        assert abs(GL.alpha_of_hwhm(GL.hwhm_of_alpha(a)) - a) < 1e-6
    a1, ak = 0.03, 0.08
    m = _microfacet_model(a1)
    out, info = GL.roughen(m, ak, kind='beckmann', shadowing=False)
    exact = _microfacet_model(np.hypot(a1, ak))
    D, De = GL.microfacet_d(out)[0], GL.microfacet_d(exact)[0]
    th = m['theta_h']; sel = np.degrees(th) < 60
    for j, k in ((0, 30), (6, 50), (12, 70)):
        assert np.nanmax(np.abs(D[j, k, sel] - De[j, k, sel])) < 1e-2 * De[j, k, 0]
        P0, P1 = (np.nan_to_num(GL._filled(X)[j, k]) * np.cos(th) ** 4 for X in (GL.microfacet_d(m)[0], D))
        r = np.tan(np.minimum(th, GL.THETA_MAX)); s = np.degrees(th) < 89
        mass = [CA.trapezoid(P[s] * 2 * np.pi * r[s], r[s]) for P in (P0, P1)]
        assert abs(mass[1] / mass[0] - 1) < 0.02
    assert np.allclose(out['S'], m['S']) and np.allclose(out['N'], m['N'])
    g = _microfacet_model(0.03, 'ggx')
    tgt, info = GL.set_gloss(g, hwhm=5.0, sym='none')
    assert abs(GL.width(tgt, 30) - 5.0) < 0.01 and abs(info['alpha_eq_after'] - GL.alpha_of_hwhm(5.0)) < 1e-6
    a1 = info['alpha_eq_after']                                                                  # finite at the horizon:
    assert np.nanmax(tgt['A']) <= np.nanmax(GL.microfacet_d(g)[0]) / a1 ** 2                     # A' = D' / (k_i k_o), k >= alpha
    noshadow, _ = GL.roughen(g, info['kernel_alpha'], shadowing=False)
    ratio = tgt['A'] / noshadow['A']
    assert np.nanmax(ratio) <= 1 + 1e-5 and abs(np.nanmean(ratio[:, 30, :20]) - 1) < 5e-3        # shadowing: only a loss, small away from grazing
    kept, info2 = GL.set_gloss(g, hwhm=5.0, sym='none', keep_albedo_deg=45)
    assert info2['albedo_after'][1] < info2['albedo_before'][1] * 1.0001 and abs(info2['albedo_after'][1] / info2['albedo_before'][1] - 1) < 1e-4
    assert info['albedo_after'][1] < info['albedo_before'][1]                                    # without it the rough lobe loses light
    try:
        GL.set_gloss(tgt, hwhm=3.0, sym='none'); raise AssertionError('a narrower target was accepted')
    except ValueError:
        pass


def test_quirks():
    """a synthetic glass with a known radiometric error K, an inpainted (copied) range below 10 deg and a one-sided arm:
    the checks give back K, the index, Fresnel's polarization with row = column, the edge of the flat range and the arm"""
    K = lambda t: 0.5 * np.cos(t) ** 2 + 0.05
    n_true = np.full(4, 1.5)
    glass = _fresnel_model(n_true, K, width_deg=0.8)
    td = glass['theta_d']; tdd = np.degrees(td); thd = np.degrees(glass['theta_h']); ph = np.degrees(glass['phi_d'])
    A = np.array(glass['A'])
    A[:, tdd < 10] = A[:, [int(np.argmin(abs(tdd - 10)))]]                            # "inpainted": a copy of the 10 deg bin
    shoulder = (thd >= 6) & (thd <= 15)
    A[:, :, shoulder] += 1e-3                                                         # a symmetric shoulder ...
    A[int(np.argmin(abs(ph - 30))), :, :][:, shoulder] *= 4                           # ... four times higher at phi_d = +30 deg
    glass['A'] = A
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'red.npz')
        QK.save(p, QK.reduce(FA.view(glass), log=None))
        m = QK.measure(QK.load(p), ('dielectric',))
    meas = (tdd >= 20) & (tdd <= 75)
    assert round(m['flat_ranges'][1e-4][0]) == 10, m['flat_ranges']                 # the copied bins change by 0
    assert abs(m['n_fit_median'] - 1.5) < 2e-3 and m['n_fit_rms'] < 2e-3
    assert np.allclose(m['energy_over_fresnel'][meas], K(td[meas]), rtol=0.02)        # the shoulder adds about 1 %
    assert np.allclose(m['N01_550'][meas], m['N01_fresnel_550'][meas], atol=2e-3) and np.allclose(m['N01_550'][meas], m['N10_550'][meas], atol=1e-6)
    a = m['arms'][25.0]
    assert round(a['phi']) == 30 and a['arm_over_mirror'] > 3.5, a
    assert abs(m['reciprocity_0_over_180'] - 1) < 1e-6
    return f"K at 45 deg {m['energy_over_fresnel'][45]:.4f} (true {K(td[45]):.4f}), arm x{a['arm_over_mirror']:.2f}"


def _body_table(n_phi=73, n_thd=31, n_thh=31, n=1.5, cos_light=1.0):
    """a diffuse body under a dielectric interface in the reader's frames: light is depolarized inside, so the Mueller matrix
    is (polarizance of the exit, by the camera's angle) x (diattenuation of the entry, by the light's angle); m00 carries
    cos(light)^cos_light, as a table fitted without the irradiance cosine does"""
    def trans(c):
        ci = np.clip(c, 1e-6, 1); ct = np.sqrt(1 - (1 - ci ** 2) / n ** 2)
        ts = 1 - ((ci - n * ct) / (ci + n * ct)) ** 2; tp = 1 - ((n * ci - ct) / (n * ci + ct)) ** 2
        return (tp - ts) / (tp + ts), 0.5 * (tp + ts)
    phi = np.linspace(-np.pi, np.pi, n_phi); td = np.linspace(0, np.pi / 2, n_thd); th = np.linspace(0, np.pi / 2, n_thh)
    TD, TH = np.meshgrid(td, th, indexing='ij')
    vals = np.full((n_phi, n_thd, n_thh, 2, 4, 4), np.nan, np.float32)
    for i, p in enumerate(phi):
        pred, cl, cv = SL.predicted(p, TD, TH)
        c1, c2, r1, r2 = pred['as labelled']
        do, ao = trans(cv); di, ai = trans(cl)
        col = np.stack([np.ones_like(c1), do * c1, do * c2, 0 * c1], -1); row = np.stack([np.ones_like(r1), di * r1, di * r2, 0 * r1], -1)
        M = (ao * ai * 0.3 * np.clip(cl, 0, 1) ** cos_light)[:, None, None] * col[:, :, None] * row[:, None, :]
        M = np.where(((cl > 0) & (cv > 0))[:, None, None], M, np.nan).reshape(TD.shape + (4, 4))
        vals[i] = M[:, :, None]
    return BRDFTable(phi, td, th, [550.0, 600.0], WAVELENGTH, True, {}, values=vals)


def test_slots():
    """Rusinkiewicz angles <-> directions; the reader's reference vectors only change sign when the directions trade places;
    a synthetic body is read 'as labelled', its exchanged view 'exchanged' (not 'mirrored'), and the one-sided cosine of m00
    follows the light"""
    rng = np.random.default_rng(0)
    a = G.normalize(rng.normal(size=(50, 3)) + [0, 0, 2]); b = G.normalize(rng.normal(size=(50, 3)) + [0, 0, 2])
    g1, g2 = G.geometry(a, b), G.geometry(b, a)
    assert np.allclose(g1['x_in'], -g2['x_out'], atol=1e-12) and np.allclose(g1['x_out'], -g2['x_in'], atol=1e-12)
    assert np.allclose(g1['theta_d'], g2['theta_d']) and np.allclose(g1['theta_h'], g2['theta_h'])
    assert np.allclose(np.abs(np.angle(np.exp(1j * (g1['phi_d'] - g2['phi_d'])))), np.pi, atol=1e-9)        # phi_d + 180 deg
    i, o = G.from_rusinkiewicz(g1['phi_d'], g1['theta_d'], g1['theta_h'])
    pd, th, td = G.rusinkiewicz(i, o)
    assert np.allclose(np.exp(1j * pd), np.exp(1j * g1['phi_d']), atol=1e-9) and np.allclose(th, g1['theta_h']) and np.allclose(td, g1['theta_d'])
    assert np.allclose(i[:, 2], a[:, 2]) and np.allclose(o[:, 2], b[:, 2])
    t = _body_table()
    kw = dict(phi_step=1, theta_h_min=0.0, theta_d_range=(0, 90), log=None)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 's.npz'); SL.save(p, SL.reduce(t, **kw)); m = SL.measure(SL.load(p), min_count=20)
    assert SL.verdict(m) == 'as labelled' and m['explained']['col']['camera'] > 0.99 and m['explained']['row']['light'] > 0.99
    assert m['alignment']['as labelled']['col'] > 0.999 and m['alignment']['exchanged']['col'] < 0.6
    assert m['cos_exponents']['light'] - m['cos_exponents']['camera'] > 0.9
    for mirrored, name, other in ((False, 'exchanged', 'exchanged, mirrored'), (True, 'exchanged, mirrored', 'exchanged')):
        mx = SL.measure(SL.reduce(CL.exchange(t, mirrored), **kw), min_count=20)
        assert SL.verdict(mx) == name and mx['explained']['col']['light'] > 0.99
        al = mx['alignment']
        assert al[name]['col'] > 0.999 and al[name]['row'] > 0.999 and al[name]['col_s2'] > 0.999, al
        assert al[other]['col_s2'] < -0.999 and al['as labelled']['col'] < 0.6, al
        assert mx['cos_exponents']['camera'] - mx['cos_exponents']['light'] > 0.9
    back = CL.exchange(CL.exchange(t))
    assert np.array_equal(np.nan_to_num(back.slice(7)), np.nan_to_num(t.slice(7)))
    return f"cos exponents light {m['cos_exponents']['light']:+.2f}, camera {m['cos_exponents']['camera']:+.2f}"


def test_clean_and_write_table():
    """dividing the light's cosine out makes the body's m00 symmetric in the two directions; g(theta_d) is divided out; a
    calibration measured with other settings is refused; a view streams into a tensor_file that reads back the same"""
    t = _body_table()
    x = CL.clean(t, exchange_slots=True, cosine='camera', min_cos=1e-3)             # exchanged: the light sits in the camera slot
    assert CL.settings(x.meta) == ('exchanged', 'camera')
    j = 9; k = int(CL.opposite_phi(t.phi_d)[j])
    a, b = x.slice(j)[..., 0, 0, 0], x.slice(k)[..., 0, 0, 0]
    ok = np.isfinite(a) & np.isfinite(b)
    assert ok.sum() > 300 and np.allclose(a[ok], b[ok], rtol=2e-4)                   # m00 now the same with light and camera exchanged
    raw_a, raw_b = t.slice(j)[..., 0, 0, 0], t.slice(k)[..., 0, 0, 0]
    assert not np.allclose(raw_a[ok], raw_b[ok], rtol=0.05)
    try:
        CL.divide_cosine(x, 'light'); raise AssertionError('a second cosine was accepted')
    except ValueError:
        pass
    g = 0.5 + 0.4 * np.cos(t.theta_d)
    y = CL.divide_theta_d(x, g)
    assert np.allclose(np.nan_to_num(y.slice(j)), np.nan_to_num(x.slice(j) / g[:, None, None, None, None].astype(np.float32)), rtol=1e-6)
    cal = dict(theta_d=t.theta_d, K=g, slots='', cosine_divided='', source='glass')
    try:
        CL.clean(t, exchange_slots=True, cosine='camera', calib=cal); raise AssertionError('a calibration with other settings was accepted')
    except ValueError:
        pass
    z = CL.clean(t, exchange_slots=True, cosine='camera', calib=dict(cal, slots='exchanged', cosine_divided='camera'), min_cos=1e-3)
    assert np.allclose(np.nan_to_num(z.slice(j)), np.nan_to_num(y.slice(j)), rtol=1e-6) and 'K(theta_d)' in CL.describe(z.meta)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'clean.pbsdf')
        hpbrdf.write_table(p, z, phi_idx=np.arange(0, t.n_phi, 2), ch_idx=[1])
        r = hpbrdf.load(p)
        assert (r.n_phi, r.n_thd, r.n_thh, r.n_ch) == (37, 31, 31, 1) and r.channels[0] == 600
        assert np.allclose(r.phi_d, t.phi_d[::2], atol=1e-6)
        assert np.array_equal(np.nan_to_num(r.slice(3)), np.nan_to_num(z.slice(6)[:, :, [1]]))
        fields = hpbrdf.parse_header(p)
        assert fields['wvls'][0] == np.uint16 and fields['phi_d'][2] == (1, 37)


def test_umov():
    """synthetic reductions: a white polarized surface + a coloured depolarized body + a scaled floor; the split gives back the
    floor's scale, the surface's Fresnel ratios, a depolarizing body and Umov's slope of -1"""
    lam = np.arange(414.0, 951.0, 8.0); td = np.radians(np.arange(0, 91, 1.0)); nT = len(td)
    c = 0.006 * (1 + 3 * ((lam - 600) / 250) ** 2)                                   # the floor's U-shaped spectrum
    B = 0.04 / (1 + np.exp(-(lam - 600) / 8))                                        # a red body: black below 560 nm
    D, C = UM.fresnel_ratios(td, 1.5)
    surf = np.zeros((nT, 4, 4)); m11 = 0.002 * (1 + np.sin(td))
    surf[:, 0, 0] = surf[:, 1, 1] = m11; surf[:, 0, 1] = surf[:, 1, 0] = D * m11; surf[:, 2, 2] = surf[:, 3, 3] = C * m11
    e00 = np.zeros((4, 4)); e00[0, 0] = 1
    body = e00.copy(); body[1, 0] = -0.05
    k_true = 0.8
    M = surf[:, None] + B[None, :, None, None] * body + (k_true * c)[None, :, None, None] * e00
    Gm = np.broadcast_to(c[None, :, None, None] * e00, M.shape) + 1e-4 * np.eye(4)   # glass: floor + a little polarized light
    n = np.full(nT, 10.0)
    red = dict(channels=lam, theta_d=td, mid_M=M * 10, mid_n=n); ref = dict(channels=lam, theta_d=td, mid_M=Gm * 10, mid_n=n)
    lam2, rows = UM.split(red, ref, dark_max=540.0)
    assert len(rows) == 11
    for r in rows:
        t = int(round(r['theta_d'])); s = r['surf']
        assert abs(r['k'] - k_true) < 1e-3 and abs(r['u']) < 1e-5, (r['k'], r['u'])
        assert abs(s[0, 1] / s[1, 1] - D[t]) < 1e-3 and abs(s[2, 2] / s[1, 1] - C[t]) < 1e-3
        assert abs(r['body'][1, 0] + 0.05) < 1e-3 and abs(r['body'][1, 1]) < 1e-3
        assert abs(r['umov_slope'] + 1) < 0.02, r['umov_slope']
    lines = []; UM.report('synthetic', lam2, rows, n=1.5, say=lines.append)
    assert len(lines) > 20


def test_predict():
    """synthetic two-component plastics: a black one (surface + floor), a white one (the same surface + a body) and a
    coloured one (body x a spectrum). The interpolation from one far-field spectrum gives the coloured one back, the tint
    does not; a target with a narrower lobe is recovered at the specular once its width is read"""
    lam = np.arange(414.0, 951.0, 8.0); td = np.radians(np.arange(0, 91, 1.0)); nT, nC = len(td), len(lam)
    c = 0.006 * (1 + 3 * ((lam - 600) / 250) ** 2)
    pig = 0.03 + 0.9 / (1 + np.exp(-(lam - 600) / 10))                                # the target's body relative to the white's
    D, C = UM.fresnel_ratios(td, 1.5)
    e00 = np.zeros((4, 4)); e00[0, 0] = 1
    nb = e00.copy(); nb[1, 0] = -0.05                                                 # the body: depolarized, a little polarizance

    def surface(level):
        S = np.zeros((nT, 4, 4)); m11 = level * (1 + np.sin(td))
        S[:, 0, 0] = S[:, 1, 1] = m11; S[:, 0, 1] = S[:, 1, 0] = D * m11; S[:, 2, 2] = S[:, 3, 3] = C * m11
        return S[:, None]

    def reduction(body, peak=1.0):
        red = dict(channels=lam, theta_d=td, theta_h=(np.arange(91) / 90.0) ** 2 * (np.pi / 2))
        for reg, s_lvl, b_lvl in (('core', 0.15 * peak, 0.05), ('mid', 0.003, 0.05), ('far', 0.0005, 0.04)):
            M = surface(s_lvl) + (b_lvl * np.cos(td))[:, None, None, None] * body[None, :, None, None] * nb + (0.8 * c)[None, :, None, None] * e00
            red[f'{reg}_M'] = M * 7.0; red[f'{reg}_n'] = np.full(nT, 7.0)
        return red
    glass = dict(channels=lam, theta_d=td)
    for reg in PR.REGIONS:
        glass[f'{reg}_M'] = np.broadcast_to(c[None, :, None, None] * e00, (nT, nC, 4, 4)) * 7.0; glass[f'{reg}_n'] = np.full(nT, 7.0)
    W, B, G = (PR.means(x) for x in (reduction(np.ones(nC)), reduction(np.zeros(nC)), glass))
    tdd = np.degrees(td)
    for peak in (1.0, 2.0):
        T = PR.means(reduction(pig, peak))
        preds, r = PR.predict(T['far'][30, :, 0, 0], W, B, G, tdd, peak=np.full(nT, peak))
        assert np.allclose(r, pig, rtol=1e-9)
        sc = PR.score(preds['interp'], T, tdd, lam)
        for reg in PR.REGIONS:
            # at the core the body's own small polarization is not modelled (2e-3 here); elsewhere the interpolation is exact
            assert sc[reg]['m00_factor'] < 1.0001 and sc[reg]['pol_rms'] < (5e-3 if reg == 'core' else 1e-6), (peak, reg, sc[reg])
        st = PR.score(preds['tint'], T, tdd, lam)
        assert st['mid']['dop_rms'] > 0.05 and st['core']['m00_rms_log'] > 0.5          # the tint dims the highlight and keeps the white's polarization
    flat, _ = PR.predict(T['far'][30, :, 0, 0], W, B, G, tdd)                                     # without the width: the peak is half too low
    assert abs(PR.score(flat['interp'], T, tdd, lam)['core']['m00_factor'] - 2.0) < 0.2
    assert abs(PR.peak_factor(np.array([30.0]), 2.0, 1.5, 4.0, 1.5)[0] - 4.0) < 1e-12
