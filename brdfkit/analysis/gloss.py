"""
Gloss as one physical knob: roughen the lobe of a factor model by convolving its distribution of microfacet slopes.

A microfacet reflection is f = F(theta_d) D(theta_h) G(i, o) / (4 cos(theta_i) cos(theta_o)). In a factor model
(m00 = A S, docs/FACTOR.md) every row of A at fixed (phi_d, theta_d) is therefore, up to the constant F(theta_d),

    D(theta_h) = A * 4 cos(theta_i) cos(theta_o)          (G ~ 1 for the polished hpBRDF materials)

and P(m) = D cos^4(theta_h), m = tan(theta_h), is the density of the microfacet slopes over the slope plane (integral 1
for a normalized D). A rougher surface of the same material has the slopes of the original plus independent slope noise,
so its density is the 2-D convolution P' = P (*) K_a, with a radially symmetric kernel: the GGX or Beckmann slope density
of roughness a. The angular part of the convolution has a closed form; the radial part is integrated on a fine grid,
folded into one [n_theta_h x n_theta_h] matrix that is applied to every row. Then

    A' = D' G'(i, o) / (4 cos(theta_i) cos(theta_o)) = D' / (k_i k_o),   k = cos + sqrt(cos^2 + a'^2 sin^2)

on the bins where A is defined, with G' the Smith shadowing of GGX at the result's equivalent roughness a' (alpha_eq): a
rough surface loses light at grazing angles, and A' stays finite at the horizon (D' / (4 cos cos) alone does not). F, S and
N are kept, and so is what the convolution only smooths, such as the haze tail of a material. The convolution preserves the
energy of the lobe; light that the wider lobe sends below the horizon is lost, as for a real rough surface in single
scattering.

The width of a lobe is the half width at half maximum (hwhm) of the phi_d-averaged D(theta_h) at a reference theta_d
(default 30 deg). It is also quoted as alpha_eq, the GGX roughness whose D has that hwhm. In the hpBRDF tables a lobe is
the instrument's response convolved with the surface's: black glass, an optically flat surface, has hwhm 1.1 deg at
theta_d 20-30 deg and about 2.5 deg at 70 deg, and every material widens with theta_d in the same way (docs/GLOSS.md).
A convolution keeps that response, so the result still looks like an hpBRDF measurement; it can only widen a lobe. A
smoother version of a material is a swap of its optics onto a sharper measured shape (`brdfkit swap`).
"""
import warnings
import numpy as np
from scipy.optimize import brentq
from scipy.special import i0e
from .factor import symmetrize

KERNELS = ('ggx', 'beckmann')
REF_THETA_D = 30.0
REPORT_THETA_D = (20.0, 30.0, 50.0, 70.0)
N_FINE = 2000                                   # fine radial nodes of the convolution (uniform in sqrt(theta_h))
THETA_MAX = np.radians(89.5)                    # slopes beyond tan(89.5 deg) are dropped
HALF_PI = np.pi / 2


# ----------------------------------------------------------------------------- GGX / Beckmann helpers
def d_ratio(theta, alpha, kind='ggx'):
    """D(theta) / D(0) of a GGX or Beckmann distribution of roughness alpha"""
    t2 = np.tan(theta) ** 2; c4 = np.cos(theta) ** 4
    if kind == 'ggx':
        return 1.0 / (c4 * (1 + t2 / alpha ** 2) ** 2)
    return np.exp(-t2 / alpha ** 2) / c4


def hwhm_of_alpha(alpha):
    """half width at half maximum [deg] of the GGX D(theta_h) of roughness alpha (alpha < 0.84)"""
    return float(np.degrees(brentq(lambda t: d_ratio(t, alpha) - 0.5, 1e-9, HALF_PI - 1e-6)))


def alpha_of_hwhm(hwhm_deg):
    """GGX roughness whose D has this half width (the inverse of hwhm_of_alpha), the lobe's alpha_eq"""
    lo, hi = 1e-5, 0.8
    if not hwhm_of_alpha(lo) < hwhm_deg < hwhm_of_alpha(hi):
        raise ValueError(f'hwhm {hwhm_deg:g} deg outside the GGX range ({hwhm_of_alpha(lo):.3g}-{hwhm_of_alpha(hi):.3g} deg)')
    return float(brentq(lambda a: hwhm_of_alpha(a) - hwhm_deg, lo, hi))


def _smith_k(c, alpha):
    """c + sqrt(c^2 + alpha^2 (1 - c^2)): the GGX Smith G1 is 2 c / k, so G1(i) G1(o) / (4 cos_i cos_o) = 1 / (k_i k_o) has no
    division by a cosine and stays finite at the horizon"""
    c = np.clip(c, 0.0, 1.0)
    return c + np.sqrt(c * c + alpha * alpha * (1.0 - c * c))


def smith_g(cos_i, cos_o, alpha):
    """separable Smith shadowing-masking of GGX, G1(i) G1(o)"""
    ci, co = np.clip(cos_i, 0.0, 1.0), np.clip(cos_o, 0.0, 1.0)
    return 4 * ci * co / (_smith_k(ci, alpha) * _smith_k(co, alpha))


# ----------------------------------------------------------------------------- geometry of the bins
def cosines(phi_d, theta_d, theta_h):
    """cos(theta_i), cos(theta_o) of every Rusinkiewicz bin [n_phi, n_thd, n_thh] (the sign convention of phi_d does not
    matter: both are used as a product or symmetrically)"""
    P, T, H = np.meshgrid(np.asarray(phi_d, float), np.asarray(theta_d, float), np.asarray(theta_h, float), indexing='ij')
    a = np.cos(H) * np.cos(T); b = np.sin(H) * np.sin(T) * np.cos(P)
    return a - b, a + b


def _applied_alpha(model):
    """alpha_eq of the Smith shadowing that `roughen` gave this model (None for a measured lobe)"""
    g = dict(model.get('meta', {})).get('gloss')
    return float(g['alpha_eq_after']) if g and g.get('shadowing') and g.get('kernel_alpha', 0) > 0 else None


def microfacet_d(model):
    """-> (D, cos_i, cos_o): D = A 4 cos_i cos_o of every bin (up to F(theta_d)), NaN where A is undefined. For a model made by
    `roughen` the shadowing it applied is divided out again (D = A k_i k_o)."""
    A = np.asarray(model['A'], np.float64)
    ci, co = cosines(model['phi_d'], model['theta_d'], model['theta_h'])
    a = _applied_alpha(model)
    D = A * 4 * ci * co if a is None else A * _smith_k(ci, a) * _smith_k(co, a)
    D[~np.isfinite(A) | (ci <= 0) | (co <= 0)] = np.nan
    return D, ci, co


def _filled(D):
    """D with its NaN bins set to the phi_d mean at the same (theta_d, theta_h), or 0: the convolution needs every slope"""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        mean = np.nan_to_num(np.nanmean(D, 0))
    return np.where(np.isfinite(D), D, mean[None])


def _hwhm(profile, theta_h):
    """hwhm [deg] of a profile over theta_h (peak at its maximum, crossing interpolated linearly)"""
    p = np.nan_to_num(np.asarray(profile, float)); th = np.degrees(np.asarray(theta_h, float))
    k = int(np.argmax(p)); half = 0.5 * p[k]
    below = np.nonzero(p[k:] < half)[0]
    if p[k] <= 0 or not len(below):
        return float('nan')
    i = k + int(below[0])
    return float(th[i - 1] + (th[i] - th[i - 1]) * (p[i - 1] - half) / (p[i - 1] - p[i]))


def _row(theta_d, deg):
    return int(np.argmin(np.abs(np.degrees(np.asarray(theta_d, float)) - deg)))


def width(model, theta_d_deg=REF_THETA_D, D=None):
    """hwhm [deg] of the phi_d-averaged D(theta_h) at the theta_d bin nearest theta_d_deg"""
    if D is None:
        D = microfacet_d(model)[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return _hwhm(np.nanmean(D[:, _row(model['theta_d'], theta_d_deg)], 0), model['theta_h'])


# ----------------------------------------------------------------------------- the convolution
def _ring(r, rp, alpha, kind):
    """integral over the angle between two slopes of length r and rp of the kernel K(|m - m'|) (closed forms)"""
    if kind == 'ggx':                                     # K(d) = a^2 / (pi (a^2 + d^2)^2)
        a = alpha ** 2 + r ** 2 + rp ** 2; b = 2 * r * rp
        return 2 * alpha ** 2 * a / ((a - b) * (a + b)) ** 1.5
    x = 2 * r * rp / alpha ** 2                           # K(d) = exp(-d^2 / a^2) / (pi a^2)
    return 2 / alpha ** 2 * i0e(x) * np.exp(-(r - rp) ** 2 / alpha ** 2)


def kernel_matrix(theta_h, alpha, kind='ggx', n_fine=N_FINE):
    """[H, H] matrix C with P' = C @ P for a slope density sampled on the theta_h nodes: P is interpolated linearly in
    sqrt(theta_h) to a fine grid, convolved, and read back at the nodes; each row is normalized to 1 (a uniform density stays
    uniform), which keeps the kernel's mass when it is narrower than the fine grid"""
    if kind not in KERNELS:
        raise ValueError(f'kernel must be one of {KERNELS}')
    th = np.asarray(theta_h, float); u = np.sqrt(np.clip(th / HALF_PI, 0, 1))
    uf = np.linspace(0.0, np.sqrt(THETA_MAX / HALF_PI), n_fine + 1)
    tf = uf ** 2 * HALF_PI; rf = np.tan(tf)
    w = np.full(len(uf), uf[1] - uf[0]); w[0] /= 2; w[-1] /= 2
    w *= rf * np.pi * uf / np.cos(tf) ** 2                 # int g(r') r' dr' with dr' = sec^2(theta) pi u du
    j = np.clip(np.searchsorted(u, uf, side='right') - 1, 0, len(u) - 2)
    f = (uf - u[j]) / (u[j + 1] - u[j])
    L = np.zeros((len(uf), len(u))); rows = np.arange(len(uf))
    L[rows, j] = 1 - f; L[rows, j + 1] = f
    r = np.tan(np.minimum(th, THETA_MAX))
    C = (_ring(r[:, None], rf[None, :], alpha, kind) * w[None, :]) @ L
    return C / C.sum(1, keepdims=True)


def convolve_d(D, theta_h, alpha, kind='ggx'):
    """D [..., H] (NaN-free) of a rougher surface: its slope density convolved with a kernel of roughness alpha"""
    th = np.asarray(theta_h, float)
    c4 = np.cos(np.minimum(th, THETA_MAX)) ** 4
    return ((D * c4) @ kernel_matrix(th, alpha, kind).T) / c4


def kernel_for(model, target_hwhm, kind='ggx', ref_theta_d=REF_THETA_D, D=None):
    """roughness of the kernel that widens the lobe to `target_hwhm` [deg] at ref_theta_d (0 if it is that wide already)"""
    if D is None:
        D = microfacet_d(model)[0]
    k = _row(model['theta_d'], ref_theta_d); th = model['theta_h']
    Dk = D[:, k]; Fk = _filled(D)[:, k]
    now = width(model, ref_theta_d, D)
    if target_hwhm <= now + 1e-3:
        if target_hwhm < now - 0.02:
            raise ValueError(f'the lobe is {now:.2f} deg wide at theta_d {ref_theta_d:g} deg; a convolution can only widen it. '
                             'For a smoother version, swap the optics onto a sharper measured shape (brdfkit swap)')
        return 0.0

    def after(a):
        Dc = convolve_d(Fk, th, a, kind)
        Dc[~np.isfinite(Dk)] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            return _hwhm(np.nanmean(Dc, 0), th)
    lo, hi = np.log(1e-4), np.log(2.0)
    if not after(np.exp(hi)) > target_hwhm:
        raise ValueError(f'target hwhm {target_hwhm:g} deg is beyond the widest kernel')
    return float(np.exp(brentq(lambda la: after(np.exp(la)) - target_hwhm, lo, hi, xtol=1e-4)))


def roughen(model, alpha_kernel, kind='ggx', shadowing=True, ref_theta_d=REF_THETA_D):
    """-> (model with its slope density convolved with a kernel of roughness alpha_kernel, info).
    D = A 4 cos_i cos_o (the shadowing of the polished originals is negligible, and dividing by it would amplify the values at
    the horizon) and, with shadowing, A' = D' G' / (4 cos_i cos_o) = D' / (k_i k_o) with G' the Smith shadowing at the
    alpha_eq after: it vanishes like cos_i cos_o at grazing angles and keeps A' finite there. Without it (G' = 1, for tests)
    A' grows without bound toward the horizon."""
    D, ci, co = microfacet_d(model)
    w0 = width(model, ref_theta_d, D)
    out = dict(model)
    if alpha_kernel > 0:
        D2 = convolve_d(_filled(D), model['theta_h'], alpha_kernel, kind)
        D2[~np.isfinite(D)] = np.nan
        w1 = width(model, ref_theta_d, D2)
        a1 = alpha_of_hwhm(w1) if shadowing else 0.0
        with np.errstate(divide='ignore', invalid='ignore'):
            A2 = D2 / (_smith_k(ci, a1) * _smith_k(co, a1))                 # D' G' / (4 cos_i cos_o)
        out['A'] = A2.astype(np.float32)
    else:
        D2 = D
    info = dict(kernel=kind, kernel_alpha=float(alpha_kernel), shadowing=bool(shadowing), ref_theta_d=float(ref_theta_d),
                hwhm_before={f'{t:g}': width(model, t, D) for t in REPORT_THETA_D},
                hwhm_after={f'{t:g}': width(model, t, D2) for t in REPORT_THETA_D})
    info.update(hwhm_ref_before=w0, hwhm_ref_after=width(model, ref_theta_d, D2))
    info.update(alpha_eq_before=alpha_of_hwhm(w0), alpha_eq_after=alpha_of_hwhm(info['hwhm_ref_after']))
    out['meta'] = dict(model.get('meta', {}), gloss=info)            # microfacet_d reads the applied shadowing from here
    return out, info


def set_gloss(model, hwhm=None, alpha=None, kernel_alpha=None, kind='ggx', shadowing=True, ref_theta_d=REF_THETA_D, sym='median',
              keep_albedo_deg=None):
    """-> (model, info) with the lobe widened to a target: hwhm [deg] or alpha (GGX alpha_eq) at ref_theta_d, or by a kernel
    of roughness kernel_alpha. A is first symmetrized over phi_d (sym, as in latent_edit.edit). The wider lobe loses the light
    that it sends below the horizon or that the shadowing removes, as a rough surface does in single scattering;
    keep_albedo_deg rescales A to keep the directional albedo at that angle of incidence instead."""
    if sum(x is not None for x in (hwhm, alpha, kernel_alpha)) != 1:
        raise ValueError('give exactly one of hwhm, alpha, kernel_alpha')
    model = symmetrize(model, sym)
    if alpha is not None:
        hwhm = hwhm_of_alpha(alpha)
    a_k = float(kernel_alpha) if kernel_alpha is not None else kernel_for(model, float(hwhm), kind, ref_theta_d)
    out, info = roughen(model, a_k, kind, shadowing, ref_theta_d)
    from .latent_edit import albedo
    scale = 1.0
    if keep_albedo_deg is not None:
        scale = albedo(model, keep_albedo_deg) / max(albedo(out, keep_albedo_deg), 1e-300)
        out['A'] = (np.asarray(out['A'], np.float64) * scale).astype(np.float32)
    info.update(target_hwhm=None if hwhm is None else float(hwhm), phi_symmetry=sym, keep_albedo_deg=keep_albedo_deg, scale=scale,
                albedo_before=[albedo(model, a) for a in (0, 45, 75)], albedo_after=[albedo(out, a) for a in (0, 45, 75)])
    out['meta'] = dict(model.get('meta', {}), gloss=info)
    return out, info


def report(info, say=print):
    ref = f'{info["ref_theta_d"]:g}'
    say(f'lobe hwhm [deg] at theta_d ' + ' / '.join(info['hwhm_before']) + ': before ' +
        ' / '.join(f'{v:.2f}' for v in info['hwhm_before'].values()) + '   after ' + ' / '.join(f'{v:.2f}' for v in info['hwhm_after'].values()))
    say(f'at theta_d {ref} deg: hwhm {info["hwhm_ref_before"]:.2f} -> {info["hwhm_ref_after"]:.2f} deg, '
        f'alpha_eq (GGX with that hwhm) {info["alpha_eq_before"]:.3f} -> {info["alpha_eq_after"]:.3f}   '
        f'({info["kernel"]} kernel, alpha {info["kernel_alpha"]:.4f}; shadowing {"on" if info["shadowing"] else "off"}; '
        f'A symmetrized over phi_d: {info.get("phi_symmetry", "none")})')
    if 'albedo_before' in info:
        say('albedo (luminance) at 0 / 45 / 75 deg: before ' + ' '.join(f'{a:.3f}' for a in info['albedo_before']) +
            '   after ' + ' '.join(f'{a:.3f}' for a in info['albedo_after']) +
            (f'   (rescaled by {info["scale"]:.3f} to keep the albedo at {info["keep_albedo_deg"]:g} deg)' if info.get('keep_albedo_deg') is not None else ''))
