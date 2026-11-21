"""
Predict a coloured plastic from a black and a white one plus a few readings of the target, and score it against its table.

A pigmented plastic reflects off its surface (white, polarized) and from inside (the pigment's colour, depolarized;
docs/UMOV.md). A black plastic is the surface alone and a white one the surface plus a full body, so at every geometry and
wavelength a plastic of another colour should lie between the two:

    M_T(geometry, lambda) = M_B + r(lambda) (M_W - M_B),      r = (R_T - R_B) / (R_W - R_B)

for all 16 Mueller elements, with r the body's weight, read from ONE spectrum R of the target at one geometry (far from the
lobe at `read_theta_d`: what a colour measurement gives). The instrument's floor is in all three tables and drops out.
The specular core needs the lobe of a sample as rough as the target, which the black one need not be. There the white
plastic's own surface is used (its m11 taken as its m00: a surface reflection keeps the polarization in the table's frames),
plus r times its body under the peak, plus the floor from black glass beside its lobe; optionally the surface is scaled by
    peak = (w_W / w_T)^2 F(n_T) / F(n_W)
from two more readings of the target, the half width w of its lobe and its index n (a lobe of the same energy and another
width, and Fresnel).

Everything works on the reductions of `brdfkit quirks`: mean Mueller matrices per theta_d and wavelength in three theta_h
regions ('core' <= 0.55 deg, 'mid' 15-25, 'far' >= 40). The baseline to beat is the tint, M_W x R_T / R_W, which is what
multiplying a white material by a colour does. docs/PREDICT.md.
"""
import numpy as np
from . import quirks as QK

REGIONS = ('core', 'mid', 'far')
ELEMENTS = ((0, 1), (1, 0), (1, 1), (2, 2), (3, 3))
K_FLOOR = 0.8                                      # the plastics' floor relative to black glass (docs/UMOV.md: 0.6-0.9)


def means(red):
    """{region: mean Mueller matrices [n_thd, n_ch, 4, 4]} of a quirks reduction"""
    out = {}
    for r in REGIONS:
        with np.errstate(invalid='ignore', divide='ignore'):
            out[r] = red[f'{r}_M'] / red[f'{r}_n'][:, None, None, None]
    return out


def fresnel_R(theta, n):
    """unpolarized Fresnel reflectance of a dielectric of index n at the angle theta [rad]"""
    ci = np.cos(theta); ct = np.sqrt(1 - (np.sin(theta) / n) ** 2)
    rs = (ci - n * ct) / (ci + n * ct); rp = (n * ci - ct) / (n * ci + ct)
    return 0.5 * (rs ** 2 + rp ** 2)


def body_weight(reading, W, B, theta_d_deg, read_theta_d=30.0):
    """r(lambda) from the target's far-field spectrum `reading` [n_ch] at read_theta_d"""
    t0 = int(np.argmin(np.abs(theta_d_deg - read_theta_d)))
    RW, RB = W['far'][t0, :, 0, 0], B['far'][t0, :, 0, 0]
    return np.clip((reading - RB) / (RW - RB), 0, None)


def peak_factor(theta_d_deg, w_target, n_target, w_white, n_white):
    th = np.radians(np.clip(theta_d_deg, 0, 89))
    return (w_white / w_target) ** 2 * fresnel_R(th, n_target) / fresnel_R(th, n_white)


def predict(reading, W, B, G, theta_d_deg, read_theta_d=30.0, peak=None, k_floor=K_FLOOR):
    """-> ({'tint': {region: M}, 'interp': {region: M}}, r). W, B, G: `means` of the white plastic, the black plastic and
    black glass; reading: the target's far-field m00 spectrum at read_theta_d; peak: None or [n_thd] (`peak_factor`)"""
    t0 = int(np.argmin(np.abs(theta_d_deg - read_theta_d)))
    r = body_weight(reading, W, B, theta_d_deg, read_theta_d)
    tint = reading / W['far'][t0, :, 0, 0]
    floor = k_floor * (G['mid'][..., 0, 0] - G['mid'][..., 1, 1])                       # under the lobe: the glass beside its own
    preds = {'tint': {}, 'interp': {}}
    for reg in REGIONS:
        preds['tint'][reg] = W[reg] * tint[None, :, None, None]
        if reg != 'core':
            preds['interp'][reg] = B[reg] + r[None, :, None, None] * (W[reg] - B[reg])
        else:
            f = np.ones(len(theta_d_deg)) if peak is None else np.asarray(peak, float)
            body = np.clip(W[reg][..., 0, 0] - floor - W[reg][..., 1, 1], 0, None)
            M = W[reg] * f[:, None, None, None]
            M[..., 0, 0] = W[reg][..., 1, 1] * f[:, None] + r[None, :] * body + floor
            preds['interp'][reg] = M
    return preds, r


def score(pred, true, theta_d_deg, channels, theta_d_pick=range(25, 80, 5), lam_max=900.0):
    """{region: dict(m00_factor, m00_rms_log, pol_rms, dop_rms)}: median |log(pred / true)| of m00 as a factor and its RMS;
    RMS error of m01, m10, m11, m22, m33 over m00; RMS error of m10 / m00 (the polarization under unpolarized light)"""
    rows = np.array([int(np.argmin(np.abs(theta_d_deg - t))) for t in theta_d_pick]); use = channels <= lam_max
    out = {}
    for reg in REGIONS:
        P = pred[reg][rows][:, use]; T = true[reg][rows][:, use]
        with np.errstate(invalid='ignore', divide='ignore'):
            lr = np.log(P[..., 0, 0] / T[..., 0, 0])
            dn = np.stack([P[..., i, j] / P[..., 0, 0] - T[..., i, j] / T[..., 0, 0] for i, j in ELEMENTS], -1)
        out[reg] = dict(m00_factor=float(np.exp(np.nanmedian(np.abs(lr)))), m00_rms_log=float(np.sqrt(np.nanmean(lr ** 2))),
                        pol_rms=float(np.sqrt(np.nanmean(dn ** 2))), dop_rms=float(np.sqrt(np.nanmean(dn[..., 1] ** 2))))
    return out


def run(target, white, black, glass, read_theta_d=30.0):
    """reductions (dicts of `quirks.load`) -> dict with the readings, r, and the scores of 'tint', 'interp' (one reading) and
    'interp+' (three readings)"""
    lam = np.asarray(target['channels'], float); td = np.degrees(np.asarray(target['theta_d'], float))
    T, W, B, G = (means(x) for x in (target, white, black, glass))
    t0 = int(np.argmin(np.abs(td - read_theta_d)))
    mT, mW = (QK.measure(x, ('dielectric',)) for x in (target, white))
    wT, wW = float(mT['hwhm_deg'][t0]), float(mW['hwhm_deg'][t0]); nT, nW = mT['n_fit_median'], mW['n_fit_median']
    peak = peak_factor(td, wT, nT, wW, nW)
    reading = T['far'][t0, :, 0, 0]
    p1, r = predict(reading, W, B, G, td, read_theta_d)
    p3, _ = predict(reading, W, B, G, td, read_theta_d, peak)
    preds = {'tint': p1['tint'], 'interp': p1['interp'], 'interp+': p3['interp']}
    return dict(channels=lam, theta_d_deg=td, r=r, readings=dict(read_theta_d=read_theta_d, hwhm=wT, n=nT, hwhm_white=wW, n_white=nW,
                peak_30=float(peak[t0]), peak_60=float(peak[int(np.argmin(np.abs(td - 60)))])),
                preds=preds, true=T, scores={k: score(p, T, td, lam) for k, p in preds.items()})


LABELS = {'tint': 'tint (white plastic x colour)', 'interp': 'interpolation, one reading (the spectrum)',
          'interp+': 'interpolation, three readings (spectrum, lobe width, index)'}


def report(name, res, say=print):
    rd = res['readings']; lam = res['channels']
    say(f'== {name}: readings of the target: far-field spectrum at theta_d {rd["read_theta_d"]:g} deg; lobe half width {rd["hwhm"]:.2f} deg '
        f'(white {rd["hwhm_white"]:.2f}); index {rd["n"]:.3f} (white {rd["n_white"]:.3f}); peak factor {rd["peak_30"]:.2f} at 30 deg, {rd["peak_60"]:.2f} at 60')
    pick = [int(np.argmin(np.abs(lam - x))) for x in (438, 502, 550, 598, 646, 702, 798, 902)]
    say('  body weight r at ' + ' '.join(f'{lam[c]:.0f}' for c in pick) + ' nm: ' + ' '.join(f'{res["r"][c]:.3f}' for c in pick))
    say('  prediction                                                    region  m00 off by (median)  RMS of log   polarization RMS   m10/m00 RMS')
    for k, sc in res['scores'].items():
        for reg in REGIONS:
            s = sc[reg]
            say(f'  {LABELS[k] if reg == REGIONS[0] else "":60s}  {reg:6s}  x{s["m00_factor"]:.2f}                {s["m00_rms_log"]:.3f}        {s["pol_rms"]:.3f}             {s["dop_rms"]:.3f}')
