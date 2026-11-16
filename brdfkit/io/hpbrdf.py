"""
Mitsuba 3 `tensor_file` tables: hpBRDF (.hpbrdf, Moon et al. 2025) and the Mitsuba export of the KAIST pBRDF (.pbsdf).

Fields: M [phi_d, theta_d, theta_h, wvl, 4, 4] float32; phi_d, theta_d, theta_h [1, n] float32 (radians); wvls [n] uint16 (nm).
Reading is lazy (one phi_d slice at a time): an hpBRDF file is 13 GB.
"""
import struct
import numpy as np
from ..table import BRDFTable, WAVELENGTH

DTYPES = {1: np.uint8, 2: np.int8, 3: np.uint16, 4: np.int16, 5: np.uint32, 6: np.int32,
          7: np.uint64, 8: np.int64, 9: np.float16, 10: np.float32, 11: np.float64}
DTYPE_CODES = {np.dtype(v): k for k, v in DTYPES.items()}
MAGIC = b'tensor_file\0'


def parse_header(fn):
    """-> {name: (dtype, byte offset, shape)}"""
    with open(fn, 'rb') as f:
        hdr = f.read(1 << 16)
    if hdr[:12] != MAGIC:
        raise ValueError(f'{fn}: not a Mitsuba tensor_file')
    p = 14
    n = struct.unpack_from('<I', hdr, p)[0]; p += 4
    fields = {}
    for _ in range(n):
        L = struct.unpack_from('<H', hdr, p)[0]; p += 2
        name = hdr[p:p + L].decode(); p += L
        nd = struct.unpack_from('<H', hdr, p)[0]; p += 2
        dt = hdr[p]; p += 1
        off = struct.unpack_from('<Q', hdr, p)[0]; p += 8
        shp = struct.unpack_from('<%dQ' % nd, hdr, p); p += 8 * nd
        fields[name] = (DTYPES[dt], off, tuple(int(s) for s in shp))
    return fields


def read_field(fn, fields, name):
    dt, off, shp = fields[name]
    return np.fromfile(fn, dtype=dt, count=int(np.prod(shp)), offset=off).reshape(shp)


def read_slice(fn, fields, i):
    """M[i] -> (n_thd, n_thh, n_wvl, 4, 4) float32"""
    dt, off, shp = fields['M']
    per = int(np.prod(shp[1:]))
    return np.fromfile(fn, dtype=dt, count=per, offset=off + np.dtype(dt).itemsize * per * i).reshape(shp[1:])


def load(fn, in_memory=False) -> BRDFTable:
    fields = parse_header(fn)
    shp = fields['M'][2]
    if len(shp) != 6 or shp[4:] != (4, 4):
        raise ValueError(f'{fn}: M has shape {shp}, expected [phi_d, theta_d, theta_h, wvl, 4, 4]')
    axes = {k: read_field(fn, fields, k).ravel().astype(np.float64) for k in ('phi_d', 'theta_d', 'theta_h')}
    wvls = read_field(fn, fields, 'wvls').ravel().astype(np.float64)
    meta = dict(source=fn, format='mitsuba tensor_file', shape=shp)
    t = BRDFTable(axes['phi_d'], axes['theta_d'], axes['theta_h'], wvls, WAVELENGTH, True, meta,
                  slice_fn=lambda i, fn=fn, fields=fields: read_slice(fn, fields, i))
    if in_memory:
        t.load_all()
    return t


def write(fn, M, phi_d, theta_d, theta_h, wvls):
    """Write a table in the same format (the fields Mitsuba's measured_polarized plugin expects)."""
    fields = {'M': np.ascontiguousarray(M, np.float32),
              'phi_d': np.ascontiguousarray(np.reshape(phi_d, (1, -1)), np.float32),
              'theta_d': np.ascontiguousarray(np.reshape(theta_d, (1, -1)), np.float32),
              'theta_h': np.ascontiguousarray(np.reshape(theta_h, (1, -1)), np.float32),
              'wvls': np.ascontiguousarray(np.ravel(wvls), np.uint16)}
    header_len = 12 + 2 + 4 + sum(2 + len(k) + 2 + 1 + 8 + 8 * v.ndim for k, v in fields.items())
    off = (header_len + 7) // 8 * 8
    offsets = {}
    for k, v in fields.items():
        offsets[k] = off
        off = (off + v.nbytes + 7) // 8 * 8
    with open(fn, 'wb') as f:
        f.write(MAGIC); f.write(struct.pack('<BB', 1, 0)); f.write(struct.pack('<I', len(fields)))
        for k, v in fields.items():
            f.write(struct.pack('<H', len(k))); f.write(k.encode())
            f.write(struct.pack('<H', v.ndim)); f.write(struct.pack('<B', DTYPE_CODES[v.dtype]))
            f.write(struct.pack('<Q', offsets[k])); f.write(struct.pack('<%dQ' % v.ndim, *v.shape))
        for k, v in fields.items():
            f.write(b'\0' * (offsets[k] - f.tell())); f.write(v.tobytes())
