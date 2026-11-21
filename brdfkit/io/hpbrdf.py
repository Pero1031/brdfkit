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


def _write_header(f, spec):
    """spec: {name: (dtype, shape)} in file order -> {name: byte offset}; leaves f at the end of the header"""
    header_len = 12 + 2 + 4 + sum(2 + len(k) + 2 + 1 + 8 + 8 * len(shp) for k, (_, shp) in spec.items())
    off = (header_len + 7) // 8 * 8
    offsets = {}
    for k, (dt, shp) in spec.items():
        offsets[k] = off
        off = (off + np.dtype(dt).itemsize * int(np.prod(shp)) + 7) // 8 * 8
    f.write(MAGIC); f.write(struct.pack('<BB', 1, 0)); f.write(struct.pack('<I', len(spec)))
    for k, (dt, shp) in spec.items():
        f.write(struct.pack('<H', len(k))); f.write(k.encode())
        f.write(struct.pack('<H', len(shp))); f.write(struct.pack('<B', DTYPE_CODES[np.dtype(dt)]))
        f.write(struct.pack('<Q', offsets[k])); f.write(struct.pack('<%dQ' % len(shp), *shp))
    return offsets


def _axes_fields(phi_d, theta_d, theta_h, wvls):
    w = np.asarray(np.ravel(wvls), float)
    if not np.allclose(w, np.round(w)) or w.min() < 0 or w.max() > 65535:
        raise ValueError('wavelengths must be whole nanometres (the format stores them as uint16)')
    return {'phi_d': np.ascontiguousarray(np.reshape(phi_d, (1, -1)), np.float32),
            'theta_d': np.ascontiguousarray(np.reshape(theta_d, (1, -1)), np.float32),
            'theta_h': np.ascontiguousarray(np.reshape(theta_h, (1, -1)), np.float32),
            'wvls': np.ascontiguousarray(np.round(w), np.uint16)}


def write(fn, M, phi_d, theta_d, theta_h, wvls):
    """Write a table in the same format (the fields Mitsuba's measured_polarized plugin expects)."""
    fields = {'M': np.ascontiguousarray(M, np.float32), **_axes_fields(phi_d, theta_d, theta_h, wvls)}
    with open(fn, 'wb') as f:
        offsets = _write_header(f, {k: (v.dtype, v.shape) for k, v in fields.items()})
        for k, v in fields.items():
            f.write(b'\0' * (offsets[k] - f.tell())); f.write(v.tobytes())


def write_table(fn, table: BRDFTable, phi_idx=None, ch_idx=None, log=None):
    """Stream a polarized spectral BRDFTable (any view, e.g. a cleaned 13 GB table or a factor model) into a tensor_file,
    one phi_d slice at a time. phi_idx / ch_idx: indices of the phi_d bins and channels to keep (default: all)."""
    if not table.polarized or table.channel_kind != WAVELENGTH:
        raise ValueError('the tensor_file format of measured_polarized needs a polarized spectral table')
    phi_idx = np.arange(table.n_phi) if phi_idx is None else np.asarray(phi_idx, int)
    ch_idx = np.arange(table.n_ch) if ch_idx is None else np.asarray(ch_idx, int)
    axes = _axes_fields(table.phi_d[phi_idx], table.theta_d, table.theta_h, table.channels[ch_idx])
    spec = {'M': (np.float32, (len(phi_idx), table.n_thd, table.n_thh, len(ch_idx), 4, 4)),
            **{k: (v.dtype, v.shape) for k, v in axes.items()}}
    all_ch = len(ch_idx) == table.n_ch and (ch_idx == np.arange(table.n_ch)).all()
    with open(fn, 'wb') as f:
        offsets = _write_header(f, spec)
        f.write(b'\0' * (offsets['M'] - f.tell()))
        for n, i in enumerate(phi_idx):
            sl = table.slice(int(i))
            f.write(np.ascontiguousarray(sl if all_ch else sl[:, :, ch_idx], '<f4').tobytes())
            if log and n % 60 == 0:
                log(f'  write {n}/{len(phi_idx)}')
        for k, v in axes.items():
            f.write(b'\0' * (offsets[k] - f.tell())); f.write(v.tobytes())
    return fn
