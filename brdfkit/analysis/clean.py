"""
Corrections of a measured table as views (each slice is corrected on access), and their composition.

What an hpBRDF table holds, as far as its own content shows (docs/QUIRKS.md, docs/SLOTS.md):

    table(phi_d, theta_d, theta_h) = [ BRDF x cos(light) x g(theta_d) ] filed at phi_d + 180 deg  +  an unpolarized floor

  floor       an additive, unpolarized, per-wavelength offset of the instrument            -> floor.subtract (docs/FLOOR.md)
  exchange    the two directions are exchanged relative to the Mitsuba reader's labels     -> exchange (docs/SLOTS.md)
  cosine      the Mueller matrices were fitted against the light as emitted, without the
              cosine of the irradiance, so the table holds BRDF x cos(light)               -> divide_cosine
  g(theta_d)  one factor per angle between light and camera (theta_d is half of it),
              measured on black glass against Fresnel (`brdfkit calibrate`)                -> divide_theta_d

`clean` applies them in that order. A calibration file must come from a glass model built with the same exchange and cosine
settings, because the K it stores is what is left after them (`brdfkit factor --exchange --cosine light`, then `calibrate`);
`settings` / `check_calibration` enforce it. Nothing here needs the 13 GB in memory.
"""
import numpy as np
from ..table import BRDFTable
from . import floor as FL

_D = np.array([1.0, 1.0, -1.0, -1.0], np.float32)
SIDES = ('light', 'camera')


def _view(table, slice_fn, **meta):
    return BRDFTable(table.phi_d, table.theta_d, table.theta_h, table.channels, table.channel_kind, table.polarized,
                     dict(table.meta, **meta), slice_fn=slice_fn)


def opposite_phi(phi_d):
    """index of the bin at phi_d + 180 deg for every bin of a full-circle phi_d axis"""
    phi = np.asarray(phi_d, float)
    wrap = (phi[None, :] - (phi[:, None] + np.pi) + np.pi) % (2 * np.pi) - np.pi
    idx = np.abs(wrap).argmin(1)
    if np.abs(wrap[np.arange(len(phi)), idx]).max() > 0.51 * np.median(np.diff(phi)):
        raise ValueError('phi_d must cover the full circle to exchange the directions')
    return idx


def exchange(table: BRDFTable, mirrored=False) -> BRDFTable:
    """view with the two directions exchanged: the slice at phi_d is the source's slice at phi_d + 180 deg. The Mueller
    matrices are kept: the reference vectors of the reader's frames only change sign when the directions trade places.
    mirrored: also flip the signs of S2 and S3 of both beams (M -> D M D, D = diag(1, 1, -1, -1))"""
    idx = opposite_phi(table.phi_d)

    def slice_fn(i):
        sl = table.slice(int(idx[i]))
        return sl * (_D[:, None] * _D[None, :]) if mirrored and table.polarized else sl
    return _view(table, slice_fn, slots='exchanged, mirrored' if mirrored else 'exchanged')


def cos_directions(table: BRDFTable, i):
    """(cos of the first direction, the reader's light; cos of the second, the camera) for the bins of phi_d slice i, [n_thd, n_thh]"""
    TD, TH = np.meshgrid(table.theta_d, table.theta_h, indexing='ij')
    a = np.cos(TD) * np.cos(TH); b = np.sin(TD) * np.cos(table.phi_d[i]) * np.sin(TH)
    return a - b, a + b


def divide_cosine(table: BRDFTable, side='light', min_cos=0.1) -> BRDFTable:
    """view with every bin divided by the cosine of its 'light' (first) or 'camera' (second) direction, clamped at min_cos:
    a table that holds BRDF x cos becomes a BRDF. The clamp keeps the last degrees before grazing (84 deg at 0.1) from
    being amplified more than 1 / min_cos"""
    if side not in SIDES:
        raise ValueError(f'side must be one of {SIDES}')
    if table.meta.get('cosine_divided'):
        raise ValueError('the table is already divided by a cosine')
    k = SIDES.index(side)

    def slice_fn(i):
        s = (1.0 / np.maximum(cos_directions(table, i)[k], min_cos)).astype(np.float32)
        sl = table.slice(i)
        return sl * s.reshape(s.shape + (1,) * (sl.ndim - 2))
    return _view(table, slice_fn, cosine_divided=side, cosine_min=float(min_cos))


def divide_theta_d(table: BRDFTable, factor, **meta) -> BRDFTable:
    """view with every bin divided by factor[theta_d]"""
    s = (1.0 / np.asarray(factor, float)).astype(np.float32)
    if s.shape != (table.n_thd,) or not np.isfinite(s).all():
        raise ValueError('factor must be finite, one value per theta_d')

    def slice_fn(i):
        sl = table.slice(i)
        return sl * s.reshape((-1,) + (1,) * (sl.ndim - 1))
    return _view(table, slice_fn, **meta)


def settings(meta):
    """(slots, cosine) of a table's or model's meta: what a calibration has to agree on"""
    return str(meta.get('slots') or ''), str(meta.get('cosine_divided') or '')


def check_calibration(meta, calib):
    have, want = settings(meta), (str(calib.get('slots', '')), str(calib.get('cosine_divided', '')))
    if have != want:
        f = lambda s: f'exchange: {"yes" if s[0] else "no"}, cosine: {s[1] or "none"}'      # noqa: E731
        raise ValueError(f'the calibration was measured with other settings ({f(want)}) than this table ({f(have)}); '
                         'build the glass model with the same --exchange / --cosine and calibrate again')


def clean(table: BRDFTable, floor=None, exchange_slots=False, cosine=None, calib=None, floor_keep=FL.KEEP,
          floor_max_gain=FL.MAX_GAIN, min_cos=0.1, mirrored=False) -> BRDFTable:
    """floor -> exchange -> cosine -> g(theta_d), each optional. floor, calib: dicts or paths; cosine: None, 'light' or 'camera'
    (labels of the table that comes out, i.e. after the exchange)"""
    if floor is not None:
        table = FL.subtract(table, floor, floor_keep, floor_max_gain)
    if exchange_slots:
        table = exchange(table, mirrored)
    if cosine:
        table = divide_cosine(table, cosine, min_cos)
    if calib is not None:
        from . import calibrate as CA
        if isinstance(calib, str):
            calib = CA.load(calib)
        check_calibration(table.meta, calib)
        table = divide_theta_d(table, CA.K_at(calib, table.theta_d), theta_d_calibration=str(calib.get('source', '')))
    return table


def describe(meta):
    """one line naming the corrections recorded in a table's meta"""
    parts = []
    if meta.get('floor_offset') is not None:
        parts.append('floor removed')
    if meta.get('slots'):
        parts.append(f'directions {meta["slots"]}')
    if meta.get('cosine_divided'):
        parts.append(f'divided by cos({meta["cosine_divided"]})')
    if 'theta_d_calibration' in meta:
        parts.append('divided by K(theta_d)')
    return ', '.join(parts) or 'as stored'
