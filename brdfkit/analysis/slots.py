"""
Which of a table's two directions is the light? A check on the diffuse body of a dielectric, and the correction it leads to.

Light that was scattered inside a dielectric leaves through its surface. Whatever went in, what comes out is polarized by the
angle it leaves at (the transmission of the interface: partly p-polarized, in the plane of the exit direction and the normal),
and how much goes in depends on the incoming polarization through the angle it arrives at. In a Mueller matrix M
(S_out = M S_in) the first is the first column (m10, m20), the polarizance, the second the first row (m01, m02). So on a
diffuse body
    |(m10, m20)| / m00  is a function of the camera's angle alone, and its orientation lies in the camera's plane of emergence;
    |(m01, m02)| / m00  is a function of the light's angle alone, and lies in the light's plane of incidence.
Neither needs a reference sample or a calibration. `reduce` bins the body of a table (far from the specular lobe) by the two
polar angles as a reader labels them and accumulates what the checks need; `measure` gives
    explained    share of the variance of each magnitude explained by the 'light' angle alone and by the 'camera' angle alone
    alignment    sum of the measured (m10, m20) projected on the predicted orientation / sum of its magnitude (+1 = p-polarized
                 as predicted, -1 = perpendicular to it), for each candidate reading of the table, and the same for the row
    cosines      exponents of  m00 = c cos^a(light) cos^b(camera)  (a BRDF's body has a = b, about 0)
The reader's labels are those of the Mitsuba `measured_polarized` plugin (docs in brdfkit/geometry.py): the first Rusinkiewicz
direction is the light. Candidate readings:
    as labelled         the first direction is the light, frames as the plugin takes them
    exchanged           the first direction is the camera; each beam keeps its reference vector x, and y = k cross x follows the
                        beam's true direction of travel. The table is then the correct one filed at phi_d + 180 deg.
    exchanged, mirrored the same with y of both beams reversed (the signs of S2 and S3 flipped)
The correction is in brdfkit/analysis/clean.py (`exchange`). docs/SLOTS.md.
"""
import json
import numpy as np
from .. import geometry as G
from ..table import BRDFTable

KIND = 'brdfkit slots reduction v1'
READINGS = ('as labelled', 'exchanged', 'exchanged, mirrored')
_N = np.array([0.0, 0.0, 1.0])


def _orientation(k_travel, toward, x):
    """unit (S1, S2) of light polarized in the plane of `toward` (the unit direction on the surface side, +-k) and the normal,
    for a beam travelling along k_travel with reference x (y = k cross x)"""
    y = np.cross(k_travel, x)
    e = _N[None, :] - toward[:, 2:3] * toward
    ex, ey = (e * x).sum(-1), (e * y).sum(-1)
    r2 = np.maximum(ex ** 2 + ey ** 2, 1e-30)
    return (ex ** 2 - ey ** 2) / r2, 2 * ex * ey / r2


def predicted(phi_d, theta_d, theta_h):
    """for bins (broadcast angles): dict reading -> (col1, col2, row1, row2), the predicted unit orientations of (m10, m20)
    and (m01, m02); plus cos of the first ('light') and second ('camera') direction"""
    i, o = G.from_rusinkiewicz(phi_d, theta_d, theta_h)
    g = G.geometry(i, o); xi, xo = g['x_in'], g['x_out']          # x of the beam arriving along -i, and of the beam leaving along +o
    out = {}
    out['as labelled'] = _orientation(o, o, xo) + _orientation(-i, i, xi)
    c1, c2 = _orientation(i, i, xi); r1, r2 = _orientation(-o, o, xo)
    out['exchanged'] = (c1, c2, r1, r2)
    out['exchanged, mirrored'] = (c1, -c2, r1, -r2)
    return out, i[:, 2], o[:, 2]


def reduce(table: BRDFTable, bands=(550.0, 700.0), phi_step=6, theta_h_min=15.0, theta_d_range=(22.0, 78.0), max_angle=80.0,
           bin_deg=10.0, log=print):
    """one pass over every `phi_step`-th phi_d slice -> dict of sums over (light angle, camera angle) bins"""
    if not table.polarized:
        raise ValueError('slots needs a polarized table')
    lam = table.channels
    bsel = (lam >= bands[0]) & (lam <= bands[1]) if table.channel_kind == 'wavelength' else np.ones(len(lam), bool)
    if not bsel.any():
        raise ValueError('no band inside --bands')
    edges = np.arange(0.0, max_angle + 1e-9, bin_deg); nb = len(edges) - 1
    tdd = np.degrees(table.theta_d)[:, None]; thd = np.degrees(table.theta_h)[None, :]
    keep = (tdd >= theta_d_range[0]) & (tdd <= theta_d_range[1]) & (thd >= theta_h_min)
    TD, TH = np.meshgrid(table.theta_d, table.theta_h, indexing='ij')
    keys = ['n', 'm00', 'col', 'row'] + [f'{w}_{k}' for w in ('col', 'row') for k in range(len(READINGS))] \
        + [f'{w}2{a}_{k}' for w in ('col', 'row') for a in ('', 'abs') for k in range(len(READINGS))]
    acc = {k: np.zeros((nb, nb)) for k in keys}
    cmin = np.cos(np.radians(max_angle))
    for i, sl in table.slices(phi_step):
        pred, cl, cv = predicted(table.phi_d[i], TD, TH)
        cl = cl.reshape(TD.shape); cv = cv.reshape(TD.shape)
        m = sl[..., bsel, :, :].astype(np.float64).mean(2)                                  # band average [T, H, 4, 4]
        ok = table.valid(sl) & keep & (cl > cmin) & (cv > cmin) & (m[..., 0, 0] > 0)
        il = np.clip(np.digitize(np.degrees(np.arccos(np.clip(cl, -1, 1))), edges) - 1, 0, nb - 1)[ok]
        iv = np.clip(np.digitize(np.degrees(np.arccos(np.clip(cv, -1, 1))), edges) - 1, 0, nb - 1)[ok]
        vals = dict(n=np.ones(ok.shape), m00=m[..., 0, 0], col=np.hypot(m[..., 1, 0], m[..., 2, 0]), row=np.hypot(m[..., 0, 1], m[..., 0, 2]))
        for k, name in enumerate(READINGS):
            c1, c2, r1, r2 = (a.reshape(TD.shape) for a in pred[name])
            vals[f'col_{k}'] = m[..., 1, 0] * c1 + m[..., 2, 0] * c2; vals[f'row_{k}'] = m[..., 0, 1] * r1 + m[..., 0, 2] * r2
            vals[f'col2_{k}'] = m[..., 2, 0] * c2; vals[f'row2_{k}'] = m[..., 0, 2] * r2
            vals[f'col2abs_{k}'] = np.abs(vals[f'col2_{k}']); vals[f'row2abs_{k}'] = np.abs(vals[f'row2_{k}'])
        for k, v in vals.items():
            np.add.at(acc[k], (il, iv), v[ok])
        if log and (i // phi_step) % 10 == 0:
            log(f'  slots {i}/{table.n_phi}')
    return dict(kind=KIND, edges=edges, **acc, meta=json.dumps(dict(
        source=str(table.meta.get('source', '')), bands=[float(lam[bsel][0]), float(lam[bsel][-1])], phi_step=phi_step,
        theta_h_min=theta_h_min, theta_d_range=list(theta_d_range), max_angle=max_angle)))


def save(path, red):
    np.savez(path, **red)


def load(path):
    with np.load(path, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    if str(d.get('kind')) != KIND:
        raise ValueError(f'{path}: not a brdfkit slots reduction')
    d['meta'] = json.loads(str(d['meta']))
    return d


def measure(red, min_count=200):
    n = red['n']; ok = n >= min_count
    if not ok.any():
        raise ValueError('no (light, camera) bin with enough samples')
    safe = np.where(red['m00'] > 0, red['m00'], 1.0)
    out = dict(ok=ok, col=np.where(ok, red['col'] / safe, np.nan), row=np.where(ok, red['row'] / safe, np.nan),
               m00=np.where(ok, red['m00'] / np.maximum(n, 1), np.nan))
    w = n * ok
    expl = {}
    for key in ('col', 'row'):
        v = np.where(ok, out[key], 0.0)
        tot = (w * (v - (w * v).sum() / w.sum()) ** 2).sum()
        e = {}
        for axis, lab in ((1, 'light'), (0, 'camera')):                                      # axis 0 of the bins is the light's angle
            mu = (w * v).sum(axis, keepdims=True) / np.maximum(w.sum(axis, keepdims=True), 1e-300)
            e[lab] = float(1 - (w * (v - mu) ** 2).sum() / max(tot, 1e-300))
        expl[key] = e
    out['explained'] = expl
    out['alignment'] = {name: dict(col=float(red[f'col_{k}'][ok].sum() / max(red['col'][ok].sum(), 1e-300)),
                                   row=float(red[f'row_{k}'][ok].sum() / max(red['row'][ok].sum(), 1e-300)),
                                   col_s2=float(red[f'col2_{k}'][ok].sum() / max(red[f'col2abs_{k}'][ok].sum(), 1e-300)),
                                   row_s2=float(red[f'row2_{k}'][ok].sum() / max(red[f'row2abs_{k}'][ok].sum(), 1e-300)))
                        for k, name in enumerate(READINGS)}
    # cosines: log m00 = log c + a log cos(light) + b log cos(camera) over the bins (centres), weighted by the counts
    e = red['edges']; cen = np.cos(np.radians(0.5 * (e[:-1] + e[1:])))
    L, V = np.meshgrid(np.log(cen), np.log(cen), indexing='ij')
    A = np.stack([np.ones(ok.sum()), L[ok], V[ok]], 1); y = np.log(out['m00'][ok]); sw = np.sqrt(w[ok])
    coef, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
    out['cos_exponents'] = dict(light=float(coef[1]), camera=float(coef[2]),
                                residual=float(np.sqrt((w[ok] * (A @ coef - y) ** 2).sum() / w[ok].sum())))
    return out


def verdict(m, margin=0.3, min_alignment=0.8):
    """one of READINGS, or None (undecided). Decided on the polarizance (the column) alone: its magnitude has to follow one
    angle and its orientation has to agree with that reading. The row is reported but not required, because a rough surface's
    reflection adds to it far from the lobe."""
    c, al = m['explained']['col'], m['alignment']
    if c['camera'] - c['light'] > margin and al['as labelled']['col'] > min_alignment:
        return 'as labelled'
    if c['light'] - c['camera'] > margin:
        for name in READINGS[1:]:
            if al[name]['col'] > min_alignment:
                return name
    return None


def report(name, red, m, say=print):
    e = red['edges']; nb = len(e) - 1; meta = red['meta']
    say(f'== {name}: body with theta_h >= {meta["theta_h_min"]:g} deg, theta_d {meta["theta_d_range"][0]:g}-{meta["theta_d_range"][1]:g} deg, '
        f'{meta["bands"][0]:.0f}-{meta["bands"][1]:.0f} nm, every {meta["phi_step"]}th phi_d slice')
    head = '   "light" \\ "camera" [deg] ' + ' '.join(f'{e[j]:3.0f}-{e[j + 1]:<3.0f}' for j in range(nb))
    for key, title, fmt, scale in (('col', 'polarizance |(m10, m20)| / m00 (the polarization of what comes out)', '{:7.3f}', 1.0),
                                   ('row', 'diattenuation |(m01, m02)| / m00 (the sensitivity to what goes in)', '{:7.3f}', 1.0),
                                   ('m00', 'm00 x 1000', '{:7.1f}', 1e3)):
        say(f'  {title}'); say(head)
        for i in range(nb):
            if m['ok'][i].any():
                say(f'   {e[i]:3.0f}-{e[i + 1]:<3.0f}                   ' + ' '.join(fmt.format(m[key][i, j] * scale) if m['ok'][i, j] else '    .  ' for j in range(nb)))
    for key, lab in (('col', 'polarizance'), ('row', 'diattenuation')):
        x = m['explained'][key]
        say(f'  {lab}: variance explained by the "light" angle alone {x["light"]:.2f}, by the "camera" angle alone {x["camera"]:.2f}')
    say('  orientation: measured vector projected on the prediction / its magnitude (+1 = p-polarized as predicted); in brackets the S2 term alone')
    for name_, a in m['alignment'].items():
        say(f'    {name_:20s} column {a["col"]:+.2f} ({a["col_s2"]:+.2f})   row {a["row"]:+.2f} ({a["row_s2"]:+.2f})')
    c = m['cos_exponents']
    say(f'  m00 = c cos^a("light") cos^b("camera"): a = {c["light"]:+.2f}, b = {c["camera"]:+.2f} (residual {c["residual"]:.2f} in log)')
    v = verdict(m)
    say('  verdict: ' + {'as labelled': 'the first direction is the light, as the reader takes it',
                         'exchanged': 'the first direction is the camera: the two directions are exchanged (the table is filed at phi_d + 180 deg)',
                         'exchanged, mirrored': 'the first direction is the camera, and S2, S3 have the other sign',
                         None: 'undecided (the polarizance does not follow one angle with the predicted orientation)'}[v])
