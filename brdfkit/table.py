"""
BRDFTable: the one in-memory representation every loader produces and every algorithm consumes.

Axes are the isotropic Rusinkiewicz angles (phi_d, theta_d, theta_h) in radians, kept exactly as the dataset samples them
(no resampling at load). Channels are either wavelengths in nm or the three RGB channels. Values are
    unpolarized : [n_phi, n_thd, n_thh, n_ch]
    polarized   : [n_phi, n_thd, n_thh, n_ch, 4, 4]   (Mueller matrices, element [0, 0] is the unpolarized BRDF)
Unmeasured bins are NaN. A table too large for memory (hpBRDF, 13 GB) provides `slice_fn(i)` returning the phi_d slice i
instead of `values`; `slice()` / `slices()` hide the difference.
"""
from dataclasses import dataclass, field
from typing import Callable, Iterator, Optional
import numpy as np

WAVELENGTH = 'wavelength'
RGB = 'rgb'


@dataclass
class BRDFTable:
    phi_d: np.ndarray
    theta_d: np.ndarray
    theta_h: np.ndarray
    channels: np.ndarray                 # wavelengths [nm] (float) or [0, 1, 2] for RGB
    channel_kind: str                    # WAVELENGTH or RGB
    polarized: bool
    meta: dict = field(default_factory=dict)
    values: Optional[np.ndarray] = None
    slice_fn: Optional[Callable[[int], np.ndarray]] = None

    def __post_init__(self):
        self.phi_d = np.asarray(self.phi_d, np.float64).ravel()
        self.theta_d = np.asarray(self.theta_d, np.float64).ravel()
        self.theta_h = np.asarray(self.theta_h, np.float64).ravel()
        self.channels = np.asarray(self.channels, np.float64).ravel()
        if self.values is None and self.slice_fn is None:
            raise ValueError('BRDFTable needs values or slice_fn')
        if self.channel_kind not in (WAVELENGTH, RGB):
            raise ValueError(f'channel_kind must be {WAVELENGTH!r} or {RGB!r}')

    # ---- sizes
    @property
    def n_phi(self): return len(self.phi_d)
    @property
    def n_thd(self): return len(self.theta_d)
    @property
    def n_thh(self): return len(self.theta_h)
    @property
    def n_ch(self): return len(self.channels)
    @property
    def slice_shape(self):
        s = (self.n_thd, self.n_thh, self.n_ch)
        return s + (4, 4) if self.polarized else s
    @property
    def in_memory(self): return self.values is not None

    # ---- access
    def slice(self, i: int) -> np.ndarray:
        """phi_d slice i -> [n_thd, n_thh, n_ch(, 4, 4)] float32"""
        if self.values is not None:
            return self.values[i]
        return self.slice_fn(i)

    def slices(self, step: int = 1) -> Iterator[tuple[int, np.ndarray]]:
        for i in range(0, self.n_phi, step):
            yield i, self.slice(i)

    def m00(self, sl: np.ndarray) -> np.ndarray:
        """unpolarized BRDF of a slice (or of any array shaped like the values)"""
        return sl[..., 0, 0] if self.polarized else sl

    def valid(self, sl: np.ndarray) -> np.ndarray:
        """mask of measured bins of a slice: finite over channels (and Mueller elements)"""
        nd = 3 if self.polarized else 1
        return np.isfinite(sl.reshape(sl.shape[:-nd] + (-1,))).all(-1)

    def load_all(self) -> np.ndarray:
        """materialize every slice (only for tables that fit in memory)"""
        if self.values is None:
            self.values = np.stack([self.slice(i) for i in range(self.n_phi)])
        return self.values

    def describe(self) -> str:
        ch = (f'{self.n_ch} wavelengths {self.channels[0]:.0f}-{self.channels[-1]:.0f} nm' if self.channel_kind == WAVELENGTH
              else 'RGB')
        return (f'BRDFTable {self.n_phi}x{self.n_thd}x{self.n_thh} bins, {ch}, '
                f'{"Mueller 4x4" if self.polarized else "unpolarized"}, {"in memory" if self.in_memory else "streamed"}')
