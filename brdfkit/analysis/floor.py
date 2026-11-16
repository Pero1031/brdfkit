"""
The instrument floor of measured tables, and a table view with it subtracted.

Far from the specular lobe the hpBRDF tables never drop below about 0.006 (550 nm), even for smooth black glass, where
the true reflectance there is close to zero. On black glass this floor
  - has the same U-shaped spectrum in all four materials checked (large where lamp x sensor sensitivity is low: 3-5x the
    550 nm value at 414 nm, 1.5-2x at 950 nm), independent of the material's own colour,
  - has no polarization structure (every m_ij / m00 below 0.02): a constant added to every raw intensity, independent of
    the QWP angles, is fitted by the Mueller reconstruction as m00 alone,
  - is flat in one of the two directions and varies about 2x with the other (not 1 / cos theta_i).
So it is read as an additive, unpolarized, per-wavelength offset of the instrument. `estimate` measures it on a reference
table (black glass); `subtract` returns a view of any table with that offset removed from m00 only.

Floor file (.npz): channels, channel_kind, offset (the statistic chosen, [n_ch]), p10, p25, p50, p90 [n_ch],
norm_absmed [n_ch, 4, 4] (median |m_ij / m00| in the region; polarized tables only), n_bins, region, source.
"""
import json
import numpy as np
from ..table import BRDFTable

KEEP = 0.02     # never take m00 below this fraction of its measured value
MAX_GAIN = 4.0  # never raise a normalized element m_ij / m00 by more than this factor (see subtract_slice)


def far_region(table: BRDFTable, theta_h_min=40.0, theta_d_range=(25.0, 75.0)):
    """[n_thd, n_thh] mask of the bins far from the specular lobe and inside the measured theta_d range of hpBRDF"""
    td = np.degrees(table.theta_d); th = np.degrees(table.theta_h)
    return ((td >= theta_d_range[0]) & (td <= theta_d_range[1]))[:, None] & (th >= theta_h_min)[None, :]


def estimate(table: BRDFTable, phi_step=12, theta_h_min=40.0, theta_d_range=(25.0, 75.0), stat='p50', log=print):
    """per-channel statistics of m00 over the far region (every `phi_step`-th phi_d slice) -> dict (see module doc)"""
    region = far_region(table, theta_h_min, theta_d_range)
    m00s, norms = [], []
    for _, sl in table.slices(phi_step):
        v = sl[region]                                                 # (n, n_ch(, 4, 4))
        flat = v.reshape(len(v), table.n_ch, -1)
        ok = np.isfinite(flat).all((1, 2)) & (table.m00(v) > 0).all(1)
        m00s.append(table.m00(v[ok]))
        if table.polarized:
            norms.append(np.abs(v[ok] / v[ok][..., :1, :1]))
    m00 = np.concatenate(m00s).astype(np.float64)
    if len(m00) == 0:
        raise ValueError('no measured bins in the far region')
    q = {f'p{p}': np.percentile(m00, p, 0) for p in (10, 25, 50, 90)}
    if stat not in q:
        raise ValueError(f'stat must be one of {sorted(q)}')
    out = dict(channels=table.channels, channel_kind=table.channel_kind, offset=q[stat], n_bins=len(m00), **q,
               region=json.dumps(dict(theta_h_min=theta_h_min, theta_d_range=list(theta_d_range), phi_step=phi_step, stat=stat)),
               source=str(table.meta.get('source', '')))
    if table.polarized:
        out['norm_absmed'] = np.median(np.concatenate(norms), 0)
    return out


def save(path, fl):
    np.savez(path, **fl)


def load(path):
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def subtract_m00(m00, offset, keep=KEEP, bound=None):
    """m00 - offset, but at least keep * m00 and at least `bound`. NaN stays NaN."""
    r = np.maximum(m00 - offset, keep * m00)
    if bound is not None:
        r = np.maximum(r, bound)
    return r


def subtract_slice(sl, offset, polarized, keep=KEEP, max_gain=MAX_GAIN):
    """one phi_d slice [n_thd, n_thh, n_ch(, 4, 4)] with the per-channel offset removed from m00.

    The offset is unpolarized, so the physically exact correction changes m00 only, which raises every m_ij / m00 by
    m00 / m00'. Where the real signal is much weaker than the offset this mostly amplifies noise; `max_gain` caps the rise:
    the other elements are scaled by min(1, max_gain * m00' / m00), so no ratio grows by more than max_gain (None = no cap)."""
    sl = np.array(sl, dtype=np.float32, copy=True)
    if not polarized:
        return subtract_m00(sl, offset, keep).astype(np.float32)
    m00 = sl[..., 0, 0].astype(np.float64)
    new = subtract_m00(m00, offset, keep)
    if max_gain is not None:
        with np.errstate(invalid='ignore', divide='ignore'):
            s = np.where(m00 > 0, np.minimum(1.0, max_gain * new / m00), 1.0)
        sl *= s[..., None, None].astype(np.float32)
    row = np.sqrt((sl[..., 0, 1:].astype(np.float64) ** 2).sum(-1)); col = np.sqrt((sl[..., 1:, 0].astype(np.float64) ** 2).sum(-1))
    absmax = np.abs(sl.reshape(sl.shape[:-2] + (16,))[..., 1:]).max(-1)                    # |m_ij| <= m00 for any Mueller matrix
    sl[..., 0, 0] = np.maximum(new, np.maximum(np.maximum(row, col), absmax))
    return sl


def subtract(table: BRDFTable, fl, keep=KEEP, max_gain=MAX_GAIN) -> BRDFTable:
    """view of `table` with the floor `fl` (dict from estimate / load, or a path) removed; slices are computed on access"""
    if isinstance(fl, str):
        fl = load(fl)
    offset = np.asarray(fl['offset'], np.float64)
    if str(fl['channel_kind']) != table.channel_kind or len(offset) != table.n_ch or \
            not np.allclose(np.asarray(fl['channels'], float), table.channels, atol=0.5):
        raise ValueError('floor channels do not match the table')
    meta = dict(table.meta, floor_offset=[float(x) for x in offset], floor_source=str(fl.get('source', '')), floor_keep=keep,
                floor_max_gain=max_gain)
    pol = table.polarized
    if table.values is not None:
        vals = np.stack([subtract_slice(table.values[i], offset, pol, keep, max_gain) for i in range(table.n_phi)])
        return BRDFTable(table.phi_d, table.theta_d, table.theta_h, table.channels, table.channel_kind, pol, meta, values=vals)
    src = table.slice_fn
    return BRDFTable(table.phi_d, table.theta_d, table.theta_h, table.channels, table.channel_kind, pol, meta,
                     slice_fn=lambda i: subtract_slice(src(i), offset, pol, keep, max_gain))


def report(fl, say=print):
    ch = np.asarray(fl['channels'], float); kind = str(fl['channel_kind'])
    reg = json.loads(str(fl['region']))
    say(f'floor of {fl["source"]}: {int(fl["n_bins"])} bins with theta_h >= {reg["theta_h_min"]} deg, '
        f'theta_d {reg["theta_d_range"][0]}-{reg["theta_d_range"][1]} deg (every {reg["phi_step"]}th phi_d slice); offset = {reg["stat"]}')
    pick = range(len(ch)) if len(ch) <= 12 else [int(np.argmin(np.abs(ch - l))) for l in (414, 450, 500, 550, 600, 650, 700, 800, 900, 950)]
    say('  channel      p10      p25      p50      p90    offset' + ('    max |m_ij/m00|' if 'norm_absmed' in fl else ''))
    for j in pick:
        lab = f'{ch[j]:.0f} nm' if kind == 'wavelength' else 'RGB'[j]
        line = f'  {lab:8s} ' + ' '.join(f'{fl[k][j]:8.5f}' for k in ('p10', 'p25', 'p50', 'p90', 'offset'))
        if 'norm_absmed' in fl:
            A = np.array(fl['norm_absmed'][j]); A[0, 0] = 0
            line += f'    {A.max():.3f}'
        say(line)
    if 'norm_absmed' in fl:
        A = np.array(fl['norm_absmed']); A[:, 0, 0] = 0
        say(f'  largest median |m_ij / m00| over all channels: {A.max():.3f} '
            '(near 0 = unpolarized offset; a reflection-like diagonal means real scattering is mixed in)')
