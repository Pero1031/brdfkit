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

DATA_HPBRDF = os.path.join(ROOT, 'data', 'kaist', 'hpbrdf', 'Aluminum.hpbrdf')
DATA_MERL = r'D:\Jack_spallow\MERL_BRDF_database\brdfs\aluminium.binary'


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


def test_merl_loader_real():
    if not os.path.exists(DATA_MERL):
        return 'skipped (no MERL file)'
    from brdfkit.io import merl
    t = merl.load(DATA_MERL)
    assert not t.polarized and t.channel_kind == RGB and t.values.shape == (180, 90, 90, 3)
    sl = t.slice(90); v = t.valid(sl)
    assert v.mean() > 0.5 and np.nanmax(sl) > 1 and np.nanmin(sl[v]) >= 0
