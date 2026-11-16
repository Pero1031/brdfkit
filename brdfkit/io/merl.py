"""
MERL BRDF database (Matusik et al. 2003): one .binary file per material.

Layout: three int32 (n_theta_h = 90, n_theta_d = 90, n_phi_d = 180), then float64[3][90 * 90 * 180], index
theta_h_idx * 90 * 180 + theta_d_idx * 180 + phi_d_idx, channel scales (1, 1.15, 1.66) / 1500. Negative values mark unmeasured bins.
Parameterization: theta_h = (i / 90)^2 * pi/2 (squared, dense near the specular peak), theta_d = i / 90 * pi/2, phi_d = i / 180 * pi
(half range: the table assumes isotropy and reciprocity). Values stored here use the bin's lower edge for each axis.
License: the 2006 package is for educational, research and non-profit use (see its Copyright_Notice.txt);
the 2023 Zenodo release (doi:10.5281/zenodo.8101680) is CC BY-SA 4.0.
"""
import numpy as np
from ..table import BRDFTable, RGB

SCALE = np.array([1.0 / 1500, 1.15 / 1500, 1.66 / 1500])


def grid_axes(n_th=90, n_td=90, n_pd=180):
    """(phi_d, theta_d, theta_h) of the MERL grid in radians (lower bin edges)"""
    return (np.arange(n_pd) / n_pd * np.pi, np.arange(n_td) / n_td * (np.pi / 2), (np.arange(n_th) / n_th) ** 2 * (np.pi / 2))


def load(fn) -> BRDFTable:
    with open(fn, 'rb') as f:
        n_th, n_td, n_pd = np.fromfile(f, np.int32, 3)
        n = int(n_th) * int(n_td) * int(n_pd)
        data = np.fromfile(f, np.float64, 3 * n)
    if data.size != 3 * n:
        raise ValueError(f'{fn}: truncated ({data.size} of {3 * n} values)')
    arr = data.reshape(3, n_th, n_td, n_pd) * SCALE[:, None, None, None]
    vals = np.ascontiguousarray(np.transpose(arr, (3, 2, 1, 0)), np.float32)       # [phi_d, theta_d, theta_h, 3]
    vals[(vals < 0).any(-1)] = np.nan
    phi_d, theta_d, theta_h = grid_axes(n_th, n_td, n_pd)
    meta = dict(source=fn, format='merl binary', license='research / educational (2006 package) or CC BY-SA 4.0 (Zenodo 2023)')
    return BRDFTable(phi_d, theta_d, theta_h, np.arange(3, dtype=np.float64), RGB, False, meta, values=vals)
