"""When groups stop being a linear superposition.

Derived in lessons/08-when-groups-go-nonlinear.md.

Everything in groups.py assumes the sea is a linear superposition of
independent components, which makes it Gaussian and its envelope Rayleigh. At
moderate steepness that stops being true: a group can pump energy into itself.
"""

import numpy as np

from .constants import G
from .dispersion import wavenumber

# Below this value of kh the sign of the nonlinear coefficient flips and
# modulational instability switches off. Benjamin & Feir's deep-water result
# simply does not apply in shallow water.
BF_CRITICAL_KH = 1.363


def nls_coefficients(omega, h=np.inf):
    """Dispersion and nonlinearity coefficients of the envelope equation

        i (A_t + c_g A_x) + alpha A_xx + beta |A|^2 A = 0.

    Deep water: alpha = -omega/(8 k^2), beta = -omega k^2 / 2.

    Returns (alpha, beta, k).
    """
    if np.isinf(h):
        k = omega**2 / G
    else:
        k = wavenumber(omega, h)
    alpha = -omega / (8.0 * k**2)
    beta = -0.5 * omega * k**2
    return alpha, beta, k


def is_focusing(k, h):
    """Focusing (unstable) above kh = 1.363, defocusing below.

    So a swell that is modulationally unstable in the open ocean becomes stable
    as it shoals. Rogue waves are a deep-water phenomenon.
    """
    return np.asarray(k) * np.asarray(h) > BF_CRITICAL_KH


def benjamin_feir_growth_rate(omega, k, a):
    """Maximum sideband growth rate, Gamma = omega (k a)^2 / 2.

    A steepness of 0.1 gives Gamma = 0.005 omega: the instability needs about
    200 wave periods to do anything, which is an hour for a 15 s swell. Slow,
    but the Pacific is wider than that.
    """
    eps = k * a
    return 0.5 * omega * eps**2


def unstable_sideband_range(k, a):
    """Perturbation wavenumbers 0 < K < 2 sqrt(2) k^2 a are unstable.

    Fastest growth at K = 2 k^2 a, i.e. a modulation wavelength of
    pi/(k^2 a) -- typically some tens of carrier wavelengths, which is exactly
    the scale of a wave group.
    """
    K_max_growth = 2.0 * k**2 * a
    K_cutoff = 2.0 * np.sqrt(2.0) * k**2 * a
    return K_max_growth, K_cutoff


def benjamin_feir_index(k_p, m0, nu):
    """BFI = sqrt(2) k_p sqrt(m0) / nu.

    The ratio of nonlinearity to spectral bandwidth. Above 1, the sea is
    steep and narrow enough for modulational instability to outrun dispersive
    spreading, and the tail of the height distribution runs above Rayleigh.

    Definitions differ. Some authors use the steepness k_p H_s / 2 rather than
    k_p sqrt(m0) = k_p H_s / 4, which shifts BFI by a factor of two. The
    threshold "about 1" should be read with that slack in mind.
    """
    if nu <= 0:
        return np.inf
    return np.sqrt(2.0) * k_p * np.sqrt(m0) / nu


def bound_infragravity_amplitude(Hs, h):
    """Order-of-magnitude amplitude of the long wave bound to a wave group.

    A group carries more radiation stress than the lull between groups, so the
    mean water level is pushed *down* beneath the big waves. That depression
    travels with the group at c_g, and when the group breaks it is set free as
    an infragravity wave with a period of 30-300 s. It is what makes the
    shorebreak surge in and out on a slow rhythm.

    This is a crude scaling, not Longuet-Higgins & Stewart's full solution:
    eta_ig ~ -(1/16) Hs^2 / h. Use it for magnitudes only.
    """
    return -(Hs**2) / (16.0 * np.asarray(h, dtype=np.float64))
