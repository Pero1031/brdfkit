"""
Directions <-> Rusinkiewicz angles, and Stokes-basis rotations.

Conventions follow the Mitsuba 3 `measured_polarized` plugin (the reference reader of the KAIST / hpBRDF tables):
  * unit vectors in the local shading frame, normal = +z
  * `to_light` points from the surface toward the light (light propagates along -to_light)
  * `to_view`  points from the surface toward the viewer (light leaves along +to_view)
  * both are rotated about z so that `to_view` has azimuth 0, then converted to (phi_d, theta_d, theta_h)
  * the table's Mueller matrix uses Stokes frames whose x reference vectors lie in the plane of reflection:
      incident beam (along -to_light): x_in;   outgoing beam (along +to_view): x_out
  * Stokes bases are right-handed (x, y = k cross x, k); the rotator sign convention is Mitsuba's / GoLD_rayt's.
All functions are vectorized over a leading sample axis: directions are (N, 3).
"""
import numpy as np

_Z = np.array([0.0, 0.0, 1.0])
_Y = np.array([0.0, 1.0, 0.0])


def normalize(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(n > 0, n, 1.0)


def rotate_z(v, ang):
    c, s = np.cos(ang), np.sin(ang)
    return np.stack([v[..., 0] * c - v[..., 1] * s, v[..., 0] * s + v[..., 1] * c, v[..., 2]], axis=-1)


def rusinkiewicz(i, o):
    """Mitsuba directions_to_rusinkiewicz(i, o) -> (phi_d, theta_h, theta_d).  i = to_light, o = to_view (std frame)."""
    h = normalize(i + o)
    b = np.cross(_Z, h)
    deg = np.linalg.norm(b, axis=-1) < 1e-7            # h == n: phi_d arbitrary
    b = np.where(deg[:, None], _Y, normalize(b))
    t = normalize(np.cross(b, h))
    td = np.arccos(np.clip((h * i).sum(-1), -1, 1))
    th = np.arccos(np.clip(h[:, 2], -1, 1))
    i_prj = normalize(i - (i * h).sum(-1, keepdims=True) * h)
    pd = np.arctan2(np.clip((b * i_prj).sum(-1), -1, 1), np.clip((t * i_prj).sum(-1), -1, 1))
    return pd, th, td


def from_rusinkiewicz(phi_d, theta_d, theta_h):
    """inverse of `rusinkiewicz` with the half vector at azimuth 0: (i, o) unit vectors, (N, 3) each (broadcast inputs)"""
    pd, td, th = np.broadcast_arrays(np.asarray(phi_d, float), np.asarray(theta_d, float), np.asarray(theta_h, float))
    pd, td, th = pd.ravel(), td.ravel(), th.ravel()
    z = np.zeros_like(th)
    h = np.stack([np.sin(th), z, np.cos(th)], -1)
    t = np.stack([np.cos(th), z, -np.sin(th)], -1)
    b = np.stack([z, z + 1.0, z], -1)
    d = np.sin(td)[:, None] * (np.cos(pd)[:, None] * t + np.sin(pd)[:, None] * b)
    return np.cos(td)[:, None] * h + d, np.cos(td)[:, None] * h - d


def geometry(to_light, to_view):
    """-> dict(phi_d, theta_d, theta_h, x_in, x_out) for (N, 3) direction arrays"""
    to_light = np.asarray(to_light, float).reshape(-1, 3); to_view = np.asarray(to_view, float).reshape(-1, 3)
    phi_std = np.arctan2(to_view[:, 1], to_view[:, 0]); phi_std = np.where(phi_std < 0, phi_std + 2 * np.pi, phi_std)
    wo = rotate_z(to_light, -phi_std); wi = rotate_z(to_view, -phi_std)
    pd, th, td = rusinkiewicz(wo, wi)
    d = wo - wi
    deg = np.linalg.norm(d, axis=-1) < 1e-7              # retro-reflection: plane of reflection undefined
    d = np.where(deg[:, None], _Y, d)
    zo = -wo; to = normalize(np.cross(d, zo)); yo = normalize(np.cross(to, zo)); xo = np.cross(yo, zo)
    zi = wi; ti = normalize(np.cross(-d, zi)); yi = normalize(np.cross(ti, zi)); xi = np.cross(yi, zi)
    return dict(phi_d=pd, theta_d=td, theta_h=th, x_in=rotate_z(xo, phi_std), x_out=rotate_z(xi, phi_std))


# ----------------------------------------------------------------------------- Stokes bases
def default_basis(w):
    """Mitsuba `stokes_basis(w)` = coordinate_system(w).first   (w: (N, 3) unit)"""
    w = np.asarray(w, float).reshape(-1, 3)
    sign = np.where(w[:, 2] >= 0, 1.0, -1.0)
    a = -1.0 / (sign + w[:, 2]); b = w[:, 0] * w[:, 1] * a
    return np.stack([w[:, 0] ** 2 * a * sign + 1.0, b * sign, -w[:, 0] * sign], axis=-1)


def rotator(theta):
    """Mueller matrix expressing the same state in a basis rotated by theta (vectorized)"""
    theta = np.asarray(theta, float)
    s, c = np.sin(2 * theta), np.cos(2 * theta)
    R = np.zeros(theta.shape + (4, 4)); R[..., 0, 0] = 1; R[..., 3, 3] = 1
    R[..., 1, 1] = c; R[..., 1, 2] = s; R[..., 2, 1] = -s; R[..., 2, 2] = c
    return R


def rotate_stokes_basis(forward, basis_current, basis_target):
    bc = normalize(basis_current); bt = normalize(basis_target)
    theta = 2 * np.arcsin(np.clip(0.5 * np.linalg.norm(bt - bc, axis=-1), -1, 1))
    theta = np.where((forward * np.cross(bc, bt)).sum(-1) < 0, -theta, theta)
    return rotator(theta)


def rotate_mueller_basis(M, in_forward, in_cur, in_tgt, out_forward, out_cur, out_tgt):
    """M expressed in (in_cur, out_cur) frames -> the same matrix in (in_tgt, out_tgt) frames"""
    R_in = rotate_stokes_basis(in_forward, in_cur, in_tgt)
    R_out = rotate_stokes_basis(out_forward, out_cur, out_tgt)
    return R_out @ M @ np.swapaxes(R_in, -1, -2)


def s_axis(to_light, to_view):
    """unit vector perpendicular to the plane of incidence (the s axis of both beams), (N, 3)"""
    s = np.cross(to_view, to_light + to_view)
    n = np.linalg.norm(s, axis=-1, keepdims=True)
    return np.where(n > 1e-9, s / np.where(n > 0, n, 1), _Y)


def mueller_to_sp(M, x_in, x_out, to_light, to_view):
    """Rotate table-frame Mueller matrices (N, 4, 4) into the (s, p) bases of both beams (s perpendicular to the plane of
    incidence, p = k cross s), the convention of Born & Wolf / GoLD_rayt. Equivalent to Mitsuba's rotate_mueller_basis."""
    s = s_axis(to_light, to_view)
    return rotate_mueller_basis(M, -to_light, x_in, s, to_view, x_out, s)
