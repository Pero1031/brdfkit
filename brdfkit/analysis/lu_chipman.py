"""
Lu-Chipman polar decomposition of Mueller matrices, Fresnel / thin-film references, and the report on a compressed pBRDF.

M = M_Delta . M_R . M_D   (Lu & Chipman, JOSA A 13, 1106 (1996))
  M_D     : diattenuator  -> diattenuation D in [0, 1]
  M_R     : retarder      -> retardance R in [0, pi]
  M_Delta : depolarizer   -> depolarization Delta in [0, 1]
plus the Gil-Bernabeu depolarization index DI (1 = non-depolarizing).

The report (`run_report`) sweeps the specular configuration over angle and wavelength, compares with the Fresnel reflection of
the bare metal (n, k from refractiveindex.info .yml files) and with one dielectric film on the metal (thickness fitted), and
maps the parameters in-plane at one wavelength. Conventions: Mueller matrices in the (s, p) bases of both beams, where an ideal
reflection at normal incidence has R = pi.
"""
import os
import numpy as np


# ----------------------------------------------------------------------------- optical constants
def load_nk(path):
    """refractiveindex.info 'tabulated nk' .yml -> (wavelength nm, n, k)"""
    rows = []
    for line in open(path, encoding='utf-8'):
        p = line.split()
        if len(p) == 3:
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                pass
    t = np.array(rows)
    return t[:, 0] * 1000.0, t[:, 1], t[:, 2]


class Metal:
    def __init__(self, path, name=None):
        self.wl, self.n, self.k = load_nk(path); self.name = name or os.path.basename(path)

    def eta(self, lam):
        return np.interp(lam, self.wl, self.n) + 1j * np.interp(lam, self.wl, self.k)


# ----------------------------------------------------------------------------- Fresnel (Born & Wolf convention)
def fresnel_rs_rp(cos_i, eta):
    """amplitudes for light arriving from index 1 onto relative index eta (complex); cos_i real in (0, 1]"""
    cos_i = np.asarray(cos_i, float); eta = np.asarray(eta, complex)
    sin2 = 1 - cos_i ** 2
    ct = np.sqrt(1 - sin2 / (eta * eta) + 0j)
    ec = eta * ct
    flip = (ec.imag < 0) | ((ec.imag == 0) & (ec.real < 0))
    ct = np.where(flip, -ct, ct); ec = eta * ct
    return (cos_i - ec) / (cos_i + ec), (eta * cos_i - ct) / (eta * cos_i + ct)


def film_rs_rp(cos_i, lam, eta_metal, n_film, d_nm):
    """one dielectric film (index n_film, thickness d) on the metal; Airy summation"""
    cos_i = np.asarray(cos_i, float)
    cos1 = np.sqrt(1 - (1 - cos_i ** 2) / n_film ** 2 + 0j)
    rs01, rp01 = fresnel_rs_rp(cos_i, n_film + 0j)
    rs12, rp12 = fresnel_rs_rp(cos1.real, eta_metal / n_film)
    e = np.exp(2j * (2 * np.pi * n_film * d_nm * cos1 / lam))
    return (rs01 + rs12 * e) / (1 + rs01 * rs12 * e), (rp01 + rp12 * e) / (1 + rp01 * rp12 * e)


def mueller_from_amplitudes(rs, rp):
    """normalized Mueller matrix (M00 = 1) of an interface in the (s, p) bases (GoLD_rayt pol::setInterface)"""
    a2, p2 = abs(rs) ** 2, abs(rp) ** 2
    nrm = 0.5 * (a2 + p2)
    c = rs * np.conj(rp) / nrm
    M = np.zeros(np.shape(rs) + (4, 4))
    M[..., 0, 0] = 1; M[..., 0, 1] = M[..., 1, 0] = 0.5 * (a2 - p2) / nrm; M[..., 1, 1] = 1
    M[..., 2, 2] = M[..., 3, 3] = c.real; M[..., 2, 3] = c.imag; M[..., 3, 2] = -c.imag
    return M


def fresnel_DR(rs, rp):
    """diattenuation and retardance of an ideal interface from its amplitudes"""
    a2, p2 = abs(rs) ** 2, abs(rp) ** 2
    return (a2 - p2) / (a2 + p2), np.abs(np.angle(rs * np.conj(rp)))


# ----------------------------------------------------------------------------- decomposition
def lu_chipman(M):
    """M (4, 4) -> dict(D, Dvec, R [rad], Raxis, Delta, Pdelta, DI, mD, mR, mDelta, detmR), or None if M00 <= 0"""
    M = np.asarray(M, float)
    m00 = M[0, 0]
    if not m00 > 0:
        return None
    M = M / m00
    Dvec = M[0, 1:].copy(); P = M[1:, 0].copy(); m = M[1:, 1:]
    Dmag = float(np.linalg.norm(Dvec))
    Dmag_c = min(Dmag, 1 - 1e-9)
    Dhat = Dvec / Dmag if Dmag > 1e-12 else np.zeros(3)
    sq = np.sqrt(1 - Dmag_c ** 2)
    mD = sq * np.eye(3) + (1 - sq) * np.outer(Dhat, Dhat)
    MD = np.eye(4); MD[0, 1:] = Dhat * Dmag_c; MD[1:, 0] = Dhat * Dmag_c; MD[1:, 1:] = mD
    mp = (M @ np.linalg.inv(MD))[1:, 1:]
    Pdelta = (P - m @ Dvec) / max(1 - Dmag_c ** 2, 1e-12)
    A = mp @ mp.T
    l1, l2, l3 = np.sqrt(np.clip(np.linalg.eigvalsh(A), 0, None))
    sign = 1.0 if np.linalg.det(mp) >= 0 else -1.0
    try:
        mDelta = sign * np.linalg.inv(A + (l1 * l2 + l2 * l3 + l3 * l1) * np.eye(3)) @ ((l1 + l2 + l3) * A + l1 * l2 * l3 * np.eye(3))
    except np.linalg.LinAlgError:
        mDelta = sign * np.eye(3) * np.trace(A) / 3
    Delta = 1 - abs(np.trace(mDelta)) / 3
    mR = np.linalg.pinv(mDelta) @ mp
    R = float(np.arccos(np.clip((np.trace(mR) - 1) / 2, -1, 1)))
    sR = np.sin(R)
    axis = (np.array([mR[1, 2] - mR[2, 1], mR[2, 0] - mR[0, 2], mR[0, 1] - mR[1, 0]]) / (2 * sR)) if abs(sR) > 1e-6 else np.full(3, np.nan)
    DI = float(np.sqrt(max((np.sum(M ** 2) - 1) / 3, 0)))
    return dict(D=Dmag, Dvec=Dvec, R=R, Raxis=axis, Delta=Delta, Pdelta=Pdelta, DI=DI, mD=mD, mR=mR, mDelta=mDelta,
                detmR=float(np.linalg.det(mR)))


def decompose_many(Ms, keys=('D', 'R', 'Delta', 'DI')):
    """(N, 4, 4) -> dict of (N,) arrays (NaN where M00 <= 0)"""
    out = {k: np.full(len(Ms), np.nan) for k in keys}
    for i, M in enumerate(Ms):
        d = lu_chipman(M)
        if d is None:
            continue
        for k in keys:
            out[k][i] = d[k]
    return out


def self_test(verbose=True):
    rng = np.random.default_rng(1)
    for th in (10, 45, 70):
        rs, rp = fresnel_rs_rp(np.cos(np.radians(th)), 1.02 + 6.63j)
        d = lu_chipman(mueller_from_amplitudes(rs, rp))
        D_ref, R_ref = fresnel_DR(rs, rp)
        assert abs(d['D'] - D_ref) < 1e-9 and abs(d['R'] - R_ref) < 1e-9 and d['Delta'] < 1e-9, (th, d, D_ref, R_ref)
    for _ in range(50):
        Dv = rng.normal(size=3); Dv *= rng.uniform(0, 0.9) / np.linalg.norm(Dv)
        Dm = np.linalg.norm(Dv); sq = np.sqrt(1 - Dm ** 2); Dh = Dv / Dm
        MD = np.eye(4); MD[0, 1:] = Dv; MD[1:, 0] = Dv; MD[1:, 1:] = sq * np.eye(3) + (1 - sq) * np.outer(Dh, Dh)
        ax = rng.normal(size=3); ax /= np.linalg.norm(ax); R = rng.uniform(0.1, 3.0)
        K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        MR = np.eye(4); MR[1:, 1:] = np.eye(3) + np.sin(R) * K + (1 - np.cos(R)) * K @ K
        diag = rng.uniform(0.3, 1.0, 3); MDl = np.eye(4); MDl[1:, 1:] = np.diag(diag)
        d = lu_chipman(0.7 * MDl @ MR @ MD)
        assert abs(d['D'] - Dm) < 1e-8 and abs(d['R'] - R) < 1e-6 and abs(d['Delta'] - (1 - diag.mean())) < 1e-8
    if verbose:
        print('lu_chipman self-test passed: Fresnel (Delta = 0, analytic D and R) and 50 synthetic products recovered')
    return True


# ----------------------------------------------------------------------------- the report
def specular_directions(thetas_deg):
    t = np.radians(np.asarray(thetas_deg, float))
    to_view = np.stack([np.sin(t), 0 * t, np.cos(t)], 1); to_light = np.stack([-np.sin(t), 0 * t, np.cos(t)], 1)
    return to_light, to_view


def principal_angle(thetas, R):
    """angle where R (rad) crosses 90 deg, linear interpolation; NaN if it never does"""
    Rd = np.degrees(R); k = np.where((Rd[:-1] >= 90) & (Rd[1:] < 90))[0]
    if len(k) == 0:
        return float('nan')
    k = k[0]; return float(thetas[k] + (Rd[k] - 90) / (Rd[k] - Rd[k + 1]) * (thetas[k + 1] - thetas[k]))


def run_report(npz, refidx, out_dir, raw=None, n_film=1.65, log=print):
    """Full analysis of a compressed polarized table. `refidx`: list of (path, name) of n,k tables (first = substrate of the film fit)."""
    from ..decode import PCADecoder
    from ..io import hpbrdf as H
    from .. import geometry as G
    os.makedirs(out_dir, exist_ok=True)
    self_test(verbose=False)
    dec = PCADecoder(npz)
    if not dec.polarized:
        raise ValueError('Lu-Chipman needs a polarized table')
    wvls = dec.channels
    metals = [Metal(p, n) for p, n in refidx]
    report = []
    def say(*a):
        s = ' '.join(str(x) for x in a); log(s); report.append(s)

    thetas = np.arange(5, 86, 2.5); t = np.radians(thetas)
    to_light, to_view = specular_directions(thetas)
    keys = ('D', 'R', 'Delta', 'DI')
    data = {k: np.full((len(thetas), len(wvls)), np.nan) for k in keys}; data['m00'] = np.full_like(data['D'], np.nan)
    for j, lam in enumerate(wvls):
        Ms = dec.eval_sp(to_light, to_view, float(lam))
        d = decompose_many(Ms, keys)
        for k in keys: data[k][:, j] = d[k]
        data['m00'][:, j] = Ms[:, 0, 0]
    fres = {}
    for met in metals:
        Dm = np.zeros_like(data['D']); Rm = np.zeros_like(data['D'])
        for j, lam in enumerate(wvls):
            Dm[:, j], Rm[:, j] = fresnel_DR(*fresnel_rs_rp(np.cos(t), met.eta(float(lam))))
        fres[met.name] = dict(D=Dm, R=Rm)

    sel = thetas <= 75; d_grid = np.arange(0, 40.01, 0.25)
    def film_err(d_nm, lam_idx):
        err = 0.0
        for j in lam_idx:
            Dm, Rm = fresnel_DR(*film_rs_rp(np.cos(t[sel]), float(wvls[j]), metals[0].eta(float(wvls[j])), n_film, d_nm))
            err += np.nanmean((Dm - data['D'][sel, j]) ** 2) + np.nanmean(((Rm - data['R'][sel, j]) / np.pi) ** 2)
        return err / len(lam_idx)
    all_idx = list(range(len(wvls)))
    errs = np.array([film_err(d, all_idx) for d in d_grid]); d_best = float(d_grid[np.argmin(errs)]); err_bare = film_err(0.0, all_idx)
    d_per_lam = [float(d_grid[np.argmin([film_err(d, [j]) for d in d_grid])]) for j in range(0, len(wvls), max(1, len(wvls) // 16))]
    film = dict(D=np.zeros_like(data['D']), R=np.zeros_like(data['D']))
    for j, lam in enumerate(wvls):
        film['D'][:, j], film['R'][:, j] = fresnel_DR(*film_rs_rp(np.cos(t), float(lam), metals[0].eta(float(lam)), n_film, d_best))

    if raw:
        try:
            fields = H.parse_header(raw)
            thd_axis = H.read_field(raw, fields, 'theta_d').ravel().astype(float)
            phi_axis = H.read_field(raw, fields, 'phi_d').ravel().astype(float)
            w550 = int(np.argmin(np.abs(wvls - 550)))
            ipi = int(np.argmin(np.abs(phi_axis - np.pi))); i0 = int(np.argmin(np.abs(phi_axis)))
            sl = H.read_slice(raw, fields, ipi)[:, 0, w550]; sl0 = H.read_slice(raw, fields, i0)[:, 0, w550]
            say(f'\nraw-table spot check at theta_h = 0, {wvls[w550]:.0f} nm (phi_d = pi slice; decoded values use the same geometry)')
            say('  theta   D raw / decoded     R raw / decoded (deg)   Delta raw / decoded   |M(phi_d=0) - M(phi_d=pi)|/m00 max')
            for th in (20, 45, 60, 75):
                jd = int(np.argmin(np.abs(thd_axis - np.radians(th)))); Mraw = sl[jd]; Mraw0 = sl0[jd]
                if not np.isfinite(Mraw).all():
                    say(f'  {th:3d}   (raw NaN)'); continue
                tl, tv = specular_directions([th]); g = G.geometry(tl, tv)
                Msp = G.mueller_to_sp(Mraw[None], g['x_in'], g['x_out'], tl, tv)[0]
                dr = lu_chipman(Msp); i = int(np.argmin(np.abs(thetas - th)))
                diff = float(np.nanmax(np.abs(Mraw0 - Mraw)) / Mraw[0, 0]) if np.isfinite(Mraw0).all() else float('nan')
                say(f'  {th:3d}   {dr["D"]:.3f} / {data["D"][i, w550]:.3f}       {np.degrees(dr["R"]):6.1f} / {np.degrees(data["R"][i, w550]):6.1f}        '
                    f'{dr["Delta"]:.3f} / {data["Delta"][i, w550]:.3f}          {diff:.3f}')
        except Exception as e:
            say('raw spot check skipped:', e)

    ang = np.arange(2.5, 86, 2.5)
    TI, TO = np.meshgrid(np.radians(ang), np.radians(ang), indexing='ij')
    tl = np.stack([-np.sin(TI), 0 * TI, np.cos(TI)], -1).reshape(-1, 3); tv = np.stack([np.sin(TO), 0 * TO, np.cos(TO)], -1).reshape(-1, 3)
    lam_map = float(wvls[int(np.argmin(np.abs(wvls - 550)))])
    Ms = dec.eval_sp(tl, tv, lam_map)
    maps = decompose_many(Ms, keys); maps['m00'] = Ms[:, 0, 0]
    maps = {k: v.reshape(TI.shape) for k, v in maps.items()}

    say(f'\n=== Lu-Chipman decomposition, specular configuration (theta_i = theta_o): {os.path.basename(npz)} vs Fresnel ===')
    say(f'film fit on {metals[0].name}: n_film = {n_film}, best thickness d = {d_best:.2f} nm (rms err {np.sqrt(errs.min()):.4f} vs bare {np.sqrt(err_bare):.4f});'
        f'  per-wavelength best d (subset of bands): {d_per_lam}')
    show_l = [l for l in (450, 550, 650, 800, 900) if wvls.min() - 1 <= l <= wvls.max() + 1] or [float(wvls[len(wvls) // 2])]
    for lam in show_l:
        j = int(np.argmin(np.abs(wvls - lam)))
        say(f'\n--- {wvls[j]:.0f} nm')
        say('  theta    D meas  D Fres  D film |  R meas  R Fres  R film (deg) |  Delta   DI     m00')
        for th in (20, 25, 30, 45, 60, 70, 80):     # hpBRDF measures theta_d ~ 20-80 deg only; below is inpainted
            i = int(np.argmin(np.abs(thetas - th)))
            say(f'  {th:3d}     {data["D"][i, j]:.3f}   {fres[metals[0].name]["D"][i, j]:.3f}   {film["D"][i, j]:.3f}  |  '
                f'{np.degrees(data["R"][i, j]):6.1f}  {np.degrees(fres[metals[0].name]["R"][i, j]):6.1f}  {np.degrees(film["R"][i, j]):6.1f}       |  '
                f'{data["Delta"][i, j]:.3f}  {data["DI"][i, j]:.3f}  {data["m00"][i, j]:.2f}')
    say('\nprincipal angle (R crosses 90 deg):')
    for lam in show_l:
        j = int(np.argmin(np.abs(wvls - lam)))
        line = f'  {wvls[j]:.0f} nm: measured {principal_angle(thetas, data["R"][:, j]):.1f} deg'
        for met in metals:
            line += f',  bare ({met.name}) {principal_angle(thetas, fres[met.name]["R"][:, j]):.1f} deg'
        say(line + f',  film {d_best:.1f} nm {principal_angle(thetas, film["R"][:, j]):.1f} deg')
    say(f'\noff-specular ({lam_map:.0f} nm, in-plane): Delta at theta_o = 10..80 step 10')
    for ti in (30, 60):
        i = int(np.argmin(np.abs(ang - ti)))
        say(f'  theta_i = {ti}:  ' + ' '.join(f'{maps["Delta"][i, int(np.argmin(np.abs(ang - to)))]:.2f}' for to in range(10, 81, 10)))

    np.savez(os.path.join(out_dir, 'lu_chipman_specular.npz'), thetas=thetas, wvls=wvls, **data,
             fresnel_D=fres[metals[0].name]['D'], fresnel_R=fres[metals[0].name]['R'], film_D=film['D'], film_R=film['R'],
             d_best=d_best, n_film=n_film, map_angles=ang, **{'map_' + k: v for k, v in maps.items()})
    open(os.path.join(out_dir, 'lu_chipman_report.txt'), 'w', encoding='utf-8').write('\n'.join(report) + '\n')
    _plots(out_dir, thetas, wvls, data, fres[metals[0].name], film, d_best, ang, maps, lam_map)
    log('wrote', out_dir)
    return dict(data=data, fresnel=fres, film=film, d_best=d_best, maps=maps)


def _plots(out_dir, thetas, wvls, data, fres, film, d_best, ang, maps, lam_map):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lams = [l for l in (450, 550, 650, 800) if wvls.min() - 1 <= l <= wvls.max() + 1] or [float(wvls[len(wvls) // 2])]
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(lams)))
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for c, lam in zip(cols, lams):
        j = int(np.argmin(np.abs(wvls - lam)))
        axs[0].plot(thetas, data['D'][:, j], '-', color=c, label=f'{wvls[j]:.0f} nm measured'); axs[0].plot(thetas, fres['D'][:, j], ':', color=c); axs[0].plot(thetas, film['D'][:, j], '--', color=c)
        axs[1].plot(thetas, np.degrees(data['R'][:, j]), '-', color=c, label=f'{wvls[j]:.0f} nm'); axs[1].plot(thetas, np.degrees(fres['R'][:, j]), ':', color=c); axs[1].plot(thetas, np.degrees(film['R'][:, j]), '--', color=c)
        axs[2].plot(thetas, data['Delta'][:, j], '-', color=c, label=f'{wvls[j]:.0f} nm Delta'); axs[2].plot(thetas, 1 - data['DI'][:, j], '-.', color=c, alpha=0.6)
    axs[0].set(title='Diattenuation D (solid: measured, dotted: bare metal, dashed: + %.1f nm film)' % d_best, xlabel='theta (deg)', ylabel='D'); axs[0].legend(fontsize=8)
    axs[1].set(title='Retardance R', xlabel='theta (deg)', ylabel='R (deg)'); axs[1].axhline(90, color='gray', lw=0.5); axs[1].legend(fontsize=8)
    axs[2].set(title='Depolarization (solid: Lu-Chipman Delta, dash-dot: 1 - DI)', xlabel='theta (deg)', ylabel='Delta'); axs[2].legend(fontsize=8)
    for a in axs: a.grid(alpha=0.3)
    fig.suptitle('specular configuration (theta_i = theta_o), (s, p) bases'); fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'lu_chipman_specular.png'), dpi=130); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    for i, th in enumerate((20, 45, 60, 75)):
        k = int(np.argmin(np.abs(thetas - th))); c = plt.cm.plasma(i / 4)
        axs[0].plot(wvls, data['D'][k], color=c, label=f'{th} deg'); axs[0].plot(wvls, fres['D'][k], ':', color=c)
        axs[1].plot(wvls, np.degrees(data['R'][k]), color=c, label=f'{th} deg'); axs[1].plot(wvls, np.degrees(fres['R'][k]), ':', color=c)
        axs[2].plot(wvls, data['Delta'][k], color=c, label=f'{th} deg')
    axs[0].set(title='D vs wavelength (dotted: bare metal)', xlabel='nm'); axs[1].set(title='R vs wavelength', xlabel='nm'); axs[2].set(title='Delta vs wavelength', xlabel='nm')
    for a in axs: a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, 'lu_chipman_spectral.png'), dpi=130); plt.close(fig)

    fig, axs = plt.subplots(1, 4, figsize=(18, 4.2)); ext = [ang[0], ang[-1], ang[0], ang[-1]]
    for a, (k, ttl, vmax) in zip(axs, (('m00', 'log10 m00', None), ('D', 'diattenuation D', 1), ('R', 'retardance R (deg)', 180), ('Delta', 'depolarization Delta', 1))):
        v = np.log10(np.clip(maps[k], 1e-3, None)) if k == 'm00' else (np.degrees(maps[k]) if k == 'R' else maps[k])
        im = a.imshow(v.T, origin='lower', extent=ext, aspect='equal', cmap='viridis' if k != 'R' else 'twilight', vmin=0 if vmax else None, vmax=vmax)
        a.plot(ang, ang, 'w--', lw=0.7); a.set(title=ttl, xlabel='theta_i (light, deg)', ylabel='theta_o (viewer, deg)'); fig.colorbar(im, ax=a, shrink=0.8)
    fig.suptitle(f'{lam_map:.0f} nm, in-plane (light and viewer on opposite sides of the normal); dashed = specular ridge')
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, f'lu_chipman_map{lam_map:.0f}.png'), dpi=130); plt.close(fig)
