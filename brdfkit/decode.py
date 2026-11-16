"""
Reference decoder for PCA coefficient tables (.npz from brdfkit.compress). The C++ decoder cpp/hpbrdf_pca.hpp implements
the same algorithm on the exported .hpbpca file; tests/ check that both agree.
"""
import json
import numpy as np
from . import geometry as G
from .compress.pca import ratio_index


def decode(c1, c2, z):
    """c1 (..., k1), c2 (..., k2) -> values (..., n_ch) or (..., n_ch, 4, 4) float32"""
    n_ch = len(z['channels'])
    m00 = np.exp(z['mean_logm00'] + c1.astype(np.float32) @ z['basis_logm00'].T)
    if not bool(z['polarized']):
        return m00.astype(np.float32)
    ratio = z['mean_ratio'] + c2.astype(np.float32) @ z['basis_ratio'].T
    out = np.empty(c1.shape[:-1] + (n_ch * 16,), np.float32)
    out[..., 0::16] = m00
    out[..., ratio_index(n_ch)] = ratio * np.repeat(m00, 15, axis=-1)
    return out.reshape(c1.shape[:-1] + (n_ch, 4, 4))


class PCADecoder:
    def __init__(self, npz_path):
        with np.load(npz_path) as zf:
            z = dict(zf)
        self.phi_d, self.theta_d, self.theta_h = (z[k].astype(np.float64) for k in ('phi_d', 'theta_d', 'theta_h'))
        self.channels = z['channels'].astype(np.float64)
        self.channel_kind = str(z['channel_kind']); self.polarized = bool(z['polarized'])
        self.k1 = z['basis_logm00'].shape[1]; self.k2 = z['basis_ratio'].shape[1]
        self.coef = np.concatenate([z['coef_logm00'], z['coef_ratio']], axis=-1).astype(np.float32)
        self.mean1 = z['mean_logm00'].astype(np.float32); self.B1 = z['basis_logm00'].astype(np.float32)
        if self.polarized:
            self.mean2 = z['mean_ratio'].reshape(-1, 15).astype(np.float32)
            self.B2 = z['basis_ratio'].reshape(-1, 15, self.k2).astype(np.float32)
        self.valid = z['valid']
        self.meta = json.loads(str(z['meta']))

    @staticmethod
    def _locate(axis, x):
        x = np.clip(x, axis[0], axis[-1])
        j = np.clip(np.searchsorted(axis, x, side='right') - 1, 0, len(axis) - 2)
        w = axis[j + 1] - axis[j]
        t = np.where(w > 0, (x - axis[j]) / np.where(w > 0, w, 1), 0.0)
        return j, t

    def coeffs(self, phi_d, theta_d, theta_h):
        """trilinear interpolation of latent coefficients -> (N, k1+k2) float32"""
        ip, tp = self._locate(self.phi_d, phi_d); id_, td = self._locate(self.theta_d, theta_d)
        ih, th = self._locate(self.theta_h, theta_h)
        c = np.zeros((len(ip), self.k1 + self.k2), np.float32)
        for dp in (0, 1):
            wp = np.where(dp, tp, 1 - tp)
            for dd in (0, 1):
                wd = np.where(dd, td, 1 - td)
                for dh in (0, 1):
                    wh = np.where(dh, th, 1 - th)
                    c += (wp * wd * wh)[:, None].astype(np.float32) * self.coef[ip + dp, id_ + dd, ih + dh]
        return c

    def _ch(self, ch):
        """channel coordinate: wavelength in nm, or 0 / 1 / 2 for RGB tables"""
        j, t = self._locate(self.channels, np.asarray(ch, float).reshape(-1))
        return j, t.astype(np.float32)

    def m00(self, c, ch):
        j, t = self._ch(ch)
        mean = (1 - t) * self.mean1[j] + t * self.mean1[j + 1]
        B = (1 - t)[:, None] * self.B1[j] + t[:, None] * self.B1[j + 1]
        return np.exp(mean + (B * c[:, :self.k1]).sum(-1))

    def mueller(self, c, ch):
        if not self.polarized:
            raise ValueError('unpolarized table')
        j, t = self._ch(ch)
        m00 = self.m00(c, ch)
        mean = (1 - t)[:, None] * self.mean2[j] + t[:, None] * self.mean2[j + 1]
        B = (1 - t)[:, None, None] * self.B2[j] + t[:, None, None] * self.B2[j + 1]
        ratio = mean + (B * c[:, None, self.k1:]).sum(-1)
        M = np.empty((len(c), 16), np.float32); M[:, 0] = m00; M[:, 1:] = ratio * m00[:, None]
        return M.reshape(-1, 4, 4)

    def spectrum_m00(self, c):
        return np.exp(self.mean1 + c[:, :self.k1] @ self.B1.T)

    def eval(self, to_light, to_view, ch, polarized=None):
        g = G.geometry(to_light, to_view)
        c = self.coeffs(g['phi_d'], g['theta_d'], g['theta_h'])
        out = dict(g, m00=self.m00(c, ch))
        if (self.polarized if polarized is None else polarized):
            out['M'] = self.mueller(c, ch)
        return out

    def eval_sp(self, to_light, to_view, ch):
        """Mueller matrices rotated into the (s, p) bases of both beams -> (N, 4, 4)"""
        to_light = np.asarray(to_light, float).reshape(-1, 3); to_view = np.asarray(to_view, float).reshape(-1, 3)
        r = self.eval(to_light, to_view, ch, polarized=True)
        return G.mueller_to_sp(r['M'], r['x_in'], r['x_out'], to_light, to_view)
