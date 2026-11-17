"""
Material swap on factor models: the angular shape of one material with the reflectance, colour and polarization of another.

A specular (microfacet) reflection is f = F(theta_d, lambda) * [D G / (4 cos_i cos_o)](angles): F holds the material's optical
constants, the bracket its surface geometry. In a factor model (docs/FACTOR.md) m00 = A S with S normalized to luminance
Y = 1, so A = level(theta_d) x shape, where the level is the luminance reflectance F_Y(theta_d) (times the instrument's
K(theta_d), the same for every hpBRDF table, docs/CALIBRATION.md), and S and N are the colour and polarization of F.
A swap keeps the shape of material X and takes from material Y, or from Fresnel of an n, k table:

  'S'      colour         S' = S_Y                       (the luminance stays: S has Y = 1)
  'N'      polarization   N' = N_Y
  'level'  reflectance    A' = A_X * level_Y / level_X   (the level is the lobe energy, calibrate.lobe_energy, smoothed in theta_d)

Any subset can be taken; all three swap the optical constants. Both models must be in the same calibration state (K cancels
in level_Y / level_X); a Fresnel source needs a calibrated X (its level is then the true reflectance). The level assumes
that the lobe holds the reflection, i.e. specular materials without a strong diffuse part (metals, glass), see docs/SWAP.md.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from . import factor as FA
from .calibrate import lobe_energy, metal_reflectance, sp_to_table
from .lu_chipman import Metal, film_rs_rp, fresnel_rs_rp, mueller_from_amplitudes

PARTS = ('S', 'N', 'level')


def _measured(model):
    td = np.degrees(np.asarray(model['theta_d'], float)); m = np.asarray(model['measured_theta_d'], float)
    return (td >= m[0] - 1e-6) & (td <= m[1] + 1e-6)


def level(model, cutoff_deg=30.0, smooth_deg=3.0):
    """luminance reflectance level of a factor model along theta_d: the lobe energy within theta_h < cutoff, smoothed over the
    measured theta_d range and copied from its edges outside"""
    q = lobe_energy(model, cutoff_deg)
    meas = _measured(model)
    if smooth_deg > 0:
        idx = np.nonzero(meas)[0]
        q = q.copy(); q[idx] = gaussian_filter1d(q[idx], smooth_deg / np.median(np.diff(np.degrees(model['theta_d']))), mode='nearest')
    return FA._extend_edges(q, meas)


def fresnel_source(model, nk_path, film_nm=0.0, n_film=1.65, name=None):
    """level, S and N of Fresnel reflection (a metal from an n, k table, optionally under one dielectric film) on the theta_d
    and wavelength grid of `model`, as a swap source"""
    td = np.asarray(model['theta_d'], float); lam = np.asarray(model['channels'], float)
    F = metal_reflectance(td, lam, nk_path, film_nm, n_film)
    met = Metal(nk_path); cos = np.cos(td)
    Msp = np.zeros((len(td), len(lam), 4, 4))
    for c, l in enumerate(lam):
        rs, rp = film_rs_rp(cos, float(l), met.eta(float(l)), n_film, film_nm) if film_nm > 0 else fresnel_rs_rp(cos, met.eta(float(l)))
        Msp[:, c] = mueller_from_amplitudes(rs, rp)
    N = sp_to_table(Msp, td).reshape(len(td), len(lam), 16)[..., 1:]
    ok = np.isfinite(N).all((1, 2)); idx = np.nonzero(ok)[0]
    N[:idx[0]] = N[idx[0]]; N[idx[-1] + 1:] = N[idx[-1]]
    yw = FA.luminance_weights(lam); FY = F @ yw
    return dict(kind='fresnel', level=FY, S=F / FY[:, None], N=N, calibrated=True,
                name=name or f'Fresnel of {nk_path}' + (f' + {film_nm:g} nm film (n {n_film:g})' if film_nm > 0 else ''))


def _calibrated(model):
    return 'radiometric_calibration' in dict(model.get('meta', {}))


def swap(shape, source, parts=PARTS, cutoff_deg=30.0, smooth_deg=3.0):
    """factor model with the angular shape of `shape` and the `parts` of `source` (a factor model or `fresnel_source`)"""
    parts = tuple(parts)
    if not parts or any(p not in PARTS for p in parts):
        raise ValueError(f'parts must be a non-empty subset of {PARTS}')
    for k in ('theta_d', 'channels'):
        if source.get('kind') != 'fresnel' and not np.allclose(shape[k], source[k]):
            raise ValueError(f'shape and source have different {k} grids')
    out = dict(shape); info = dict(parts=list(parts), cutoff_deg=cutoff_deg, smooth_deg=smooth_deg)
    if 'level' in parts:
        if source.get('kind') == 'fresnel':
            if not _calibrated(shape):
                raise ValueError('a Fresnel source needs a calibrated shape model (brdfkit edit --calib)')
            lev_src = np.asarray(source['level'], float)
        else:
            if _calibrated(shape) != _calibrated(source):
                raise ValueError('calibrate both models or neither: the instrument error K cancels only between equals')
            lev_src = level(source, cutoff_deg, smooth_deg)
        lev_shape = level(shape, cutoff_deg, smooth_deg)
        ratio = lev_src / np.maximum(lev_shape, 1e-12)
        out['A'] = (np.asarray(shape['A'], np.float32) * ratio[None, :, None].astype(np.float32))
        info['level_ratio_median'] = float(np.median(ratio[_measured(shape)]))
    if 'S' in parts:
        out['S'] = np.asarray(source['S'], float).copy()
    if 'N' in parts:
        out['N'] = np.asarray(source['N'], float).copy()
    src_name = source.get('name') if source.get('kind') == 'fresnel' else str(dict(source.get('meta', {})).get('source', ''))
    info.update(shape=str(dict(shape.get('meta', {})).get('source', '')), source=src_name)
    out['meta'] = dict(shape.get('meta', {}), swap=info)
    return out, info


def report(info, say=print):
    say(f'shape of {info["shape"]}; {", ".join(info["parts"])} from {info["source"]}' +
        (f'; reflectance level x {info["level_ratio_median"]:.3g} (median over theta_d)' if 'level_ratio_median' in info else ''))
