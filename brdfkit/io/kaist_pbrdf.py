"""
KAIST pBRDF dataset (Baek et al., SIGGRAPH 2020): MATLAB .mat with pbrdf_T [phi_d 361, theta_d 91, theta_h 91, wvl 5, 4, 4],
pbrdf_T_f (inpainted), phi_d / theta_d / theta_h [1, n] (radians), wvls [1, 5] uint16 (nm). About 1.9 GB per material in memory.
"""
import numpy as np
from ..table import BRDFTable, WAVELENGTH


def load(fn, inpainted=False) -> BRDFTable:
    import scipy.io as sio
    key = 'pbrdf_T_f' if inpainted else 'pbrdf_T'
    d = sio.loadmat(fn, variable_names=[key, 'phi_d', 'theta_d', 'theta_h', 'wvls'])
    M = np.ascontiguousarray(d[key], np.float32)
    meta = dict(source=fn, format='kaist pbrdf .mat', variable=key,
                license='no formal license; the dataset page asks for a citation of Baek et al. 2020')
    return BRDFTable(d['phi_d'].ravel(), d['theta_d'].ravel(), d['theta_h'].ravel(), d['wvls'].ravel().astype(np.float64),
                     WAVELENGTH, True, meta, values=M)
