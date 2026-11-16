"""Dataset loaders. `load(path)` picks one by file extension (.npz: a factor model of `brdfkit factor`, rebuilt on access)."""
import os
from ..table import BRDFTable


def load(path: str, **kw) -> BRDFTable:
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.hpbrdf', '.pbsdf'):
        from . import hpbrdf
        return hpbrdf.load(path, **kw)
    if ext == '.mat':
        from . import kaist_pbrdf
        return kaist_pbrdf.load(path, **kw)
    if ext == '.binary':
        from . import merl
        return merl.load(path, **kw)
    if ext == '.npz':
        from ..analysis import factor
        if factor.is_model(path):
            return factor.view(path)
    raise ValueError(f'unknown BRDF file type: {path}')
