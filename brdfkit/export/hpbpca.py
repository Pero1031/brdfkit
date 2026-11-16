"""
.hpbpca v2: the flat binary the C++ decoder (cpp/hpbrdf_pca.hpp) reads. Little endian.

  char    magic[8] = "HPBPCA02"
  uint32  n_phi, n_thd, n_thh, n_ch, k1, k2, coef_type (0 = float16, 1 = float32), flags
          flags bit 0: polarized (Mueller part present);  bit 1: channels are RGB (else wavelengths in nm)
  float   eps
  float   phi_d[n_phi], theta_d[n_thd], theta_h[n_thh]            radians, ascending
  float   channels[n_ch]                                          nm, or 0 1 2 for RGB
  float   mean_logm00[n_ch], basis_logm00[n_ch][k1]
  float   mean_ratio[n_ch*15], basis_ratio[n_ch*15][k2]           only when polarized (ratio index = ch*15 + (e-1), e = row*4+col)
  uint8   valid[n_phi][n_thd][n_thh]
  coef_t  coef[n_phi][n_thd][n_thh][k1+k2]
"""
import json
import os
import struct
import numpy as np
from ..table import RGB

MAGIC = b'HPBPCA02'


def write(npz_path, out, dtype='float16'):
    with np.load(npz_path) as zf:
        z = dict(zf)
    C1, C2 = z['coef_logm00'], z['coef_ratio']
    n_phi, n_thd, n_thh, k1 = C1.shape; k2 = C2.shape[-1]; n_ch = len(z['channels'])
    polarized = bool(z['polarized']); rgb = str(z['channel_kind']) == RGB
    dt = np.float32 if dtype == 'float32' else np.float16
    coef = np.concatenate([C1, C2], axis=-1).astype(dt)
    flags = (1 if polarized else 0) | (2 if rgb else 0)
    with open(out, 'wb') as f:
        f.write(MAGIC)
        f.write(struct.pack('<8I', n_phi, n_thd, n_thh, n_ch, k1, k2, 1 if dt == np.float32 else 0, flags))
        f.write(struct.pack('<f', float(json.loads(str(z['meta']))['eps'])))
        for k in ('phi_d', 'theta_d', 'theta_h', 'channels'):
            f.write(np.ascontiguousarray(z[k], np.float32).tobytes())
        f.write(np.ascontiguousarray(z['mean_logm00'], np.float32).tobytes())
        f.write(np.ascontiguousarray(z['basis_logm00'], np.float32).tobytes())
        if polarized:
            f.write(np.ascontiguousarray(z['mean_ratio'], np.float32).tobytes())
            f.write(np.ascontiguousarray(z['basis_ratio'], np.float32).tobytes())
        f.write(np.ascontiguousarray(z['valid'], np.uint8).tobytes())
        f.write(coef.tobytes())
    size = os.path.getsize(out)
    print(f'wrote {out}  {size/2**20:.1f} MB  k1={k1} k2={k2} coef={dt.__name__} polarized={polarized} rgb={rgb}')
    return out


def read(path):
    """-> dict with the same keys as the compressor's .npz (coefficients split back into coef_logm00 / coef_ratio)"""
    with open(path, 'rb') as f:
        if f.read(8) != MAGIC:
            raise ValueError(f'{path}: not an HPBPCA02 file')
        n_phi, n_thd, n_thh, n_ch, k1, k2, coef_type, flags = struct.unpack('<8I', f.read(32))
        eps = struct.unpack('<f', f.read(4))[0]
        rd = lambda n, dt=np.float32: np.fromfile(f, dt, n)
        d = dict(phi_d=rd(n_phi), theta_d=rd(n_thd), theta_h=rd(n_thh), channels=rd(n_ch))
        d['mean_logm00'] = rd(n_ch); d['basis_logm00'] = rd(n_ch * k1).reshape(n_ch, k1)
        polarized = bool(flags & 1)
        if polarized:
            d['mean_ratio'] = rd(n_ch * 15); d['basis_ratio'] = rd(n_ch * 15 * k2).reshape(n_ch * 15, k2)
        else:
            d['mean_ratio'] = np.zeros(0, np.float32); d['basis_ratio'] = np.zeros((0, 0), np.float32)
        d['valid'] = rd(n_phi * n_thd * n_thh, np.uint8).reshape(n_phi, n_thd, n_thh).astype(bool)
        coef = rd((n_phi * n_thd * n_thh) * (k1 + k2), np.float32 if coef_type == 1 else np.float16).reshape(n_phi, n_thd, n_thh, k1 + k2)
    d['coef_logm00'] = coef[..., :k1]; d['coef_ratio'] = coef[..., k1:]
    d['polarized'] = polarized; d['channel_kind'] = RGB if flags & 2 else 'wavelength'
    d['meta'] = json.dumps(dict(eps=eps, k1=k1, k2=k2, source=path))
    return d
