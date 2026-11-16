"""
Spectral band weights: wavelength samples -> linear sRGB.

Identical to the RGB bake of `integrations/gold_rayt/include/Materials/MeasuredHpBRDF.hpp` (a test compares the constants):
CIE 1931 2-degree colour matching functions as the multi-lobe Gaussian fit of Wyman, Sloan, Shirley (2013), XYZ -> linear sRGB
(D65 matrix), each channel's weights normalized to sum to 1 over the bands so that a flat spectrum maps to equal RGB.
The weights can be slightly negative (sRGB gamut); callers clamp the resulting RGB at 0 as the renderer does.
"""
import numpy as np

# (amplitude, centre [nm], sigma below the centre, sigma above the centre) of each lobe
CMF_LOBES = {
    'x': ((1.056, 599.8, 37.9, 31.0), (0.362, 442.0, 16.0, 26.7), (-0.065, 501.1, 20.4, 26.2)),
    'y': ((0.821, 568.8, 46.9, 40.5), (0.286, 530.9, 16.3, 31.1)),
    'z': ((1.217, 437.0, 11.8, 36.0), (0.681, 459.0, 26.0, 13.8)),
}
XYZ_TO_SRGB = np.array([[3.2406, -1.5372, -0.4986],
                        [-0.9689, 1.8758, 0.0415],
                        [0.0557, -0.2040, 1.0570]])


def cmf(lam):
    """wavelengths [nm] (n,) -> (n, 3) CIE XYZ matching functions"""
    lam = np.asarray(lam, np.float64)
    out = []
    for axis in 'xyz':
        s = np.zeros_like(lam)
        for a, mu, s1, s2 in CMF_LOBES[axis]:
            t = (lam - mu) / np.where(lam < mu, s1, s2)
            s += a * np.exp(-0.5 * t * t)
        out.append(s)
    return np.stack(out, -1)


def rgb_band_weights(wavelengths):
    """-> (n_bands, 3) weights: rgb = spectrum @ weights, each column sums to 1"""
    w = cmf(wavelengths) @ XYZ_TO_SRGB.T
    return w / w.sum(0, keepdims=True)
