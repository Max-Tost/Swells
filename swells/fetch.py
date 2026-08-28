"""Fetch- and duration-limited wave growth.

Derived in lessons/04-the-spectrum.md.

Everything in this module is empirical. The scaling laws come from the JONSWAP
field campaign in the North Sea (Hasselmann et al. 1973) and the duration
relation from the CERC Shore Protection Manual lineage (Carter 1982). They are
fits to data, not consequences of the equations of motion, and they are quoted
here with the ranges over which they were fitted.
"""

from dataclasses import dataclass

import numpy as np

from .constants import G

# Dimensionless fetch at which the sea stops growing and becomes "fully
# developed". Obtained by setting the JONSWAP peak-frequency law equal to the
# Pierson-Moskowitz value f_p U10 / g = 0.13:
#     3.5 * x^-0.33 = 0.13  =>  x = (3.5/0.13)^(1/0.33) = 2.16e4
FULLY_DEVELOPED_FETCH = 2.161e4

# Validity range of the JONSWAP fit (Hasselmann et al. 1973, Fig. 2.3).
JONSWAP_FETCH_RANGE = (1e2, 1e5)

# JONSWAP total-energy growth law:  eps = g^2 m0 / U10^4 = 1.6e-7 x~
# Equivalently g H_s / U10^2 = 0.0016 sqrt(x~), since H_s = 4 sqrt(m0).
ENERGY_GROWTH_COEFF = 1.6e-7

# THE INCONSISTENCY, stated plainly because it bit the draft this project was
# built from.
#
# JONSWAP gives three separate fits: for alpha, for f_p, and for total energy.
# They do not agree with each other. Integrating the spectrum built from the
# alpha and f_p fits gives
#
#     g H_s / U^2 = 1.26e-3 x~^0.55
#
# while the measured energy fit says
#
#     g H_s / U^2 = 1.60e-3 x~^0.50
#
# Different coefficient and different exponent. At x~ = 1e4 the spectral route
# comes out 24% high; at x~ = 1e3 it is 4% low. This is not a bug in anyone's
# algebra -- they are three independent regressions through scattered field
# data, and nobody made them close.
#
# We resolve it by trusting the energy law, because m_0 is the thing that was
# measured most directly, and rescaling the spectrum to match. The shape stays
# JONSWAP; only the level moves. Set calibrate=False to get the raw spectrum
# and see the discrepancy for yourself (tests/test_spectra.py measures it).


@dataclass
class SeaState:
    """What the wind has managed to build."""
    U10: float                 # wind speed at 10 m, m/s
    fetch: float               # geometric fetch, m
    duration: float            # wind duration, s (inf if unconstrained)
    x_tilde: float             # dimensionless fetch actually governing growth
    f_p: float                 # peak frequency, Hz
    alpha: float               # Phillips-type scale constant
    regime: str                # 'fetch-limited', 'duration-limited', 'fully-developed'
    energy_scale: float = 1.0  # factor applied to reconcile the two fits

    @property
    def T_p(self):
        return 1.0 / self.f_p

    @property
    def target_m0(self):
        """m_0 demanded by the JONSWAP energy law."""
        return ENERGY_GROWTH_COEFF * self.x_tilde * self.U10**4 / G**2

    @property
    def Hs(self):
        return 4.0 * np.sqrt(self.target_m0)


def dimensionless_fetch(fetch, U10):
    """x~ = g F / U10^2."""
    return G * fetch / U10**2


def dimensionless_duration(duration, U10):
    """t~ = g t / U10."""
    return G * duration / U10


def duration_limited_fetch(duration, U10):
    """The fetch a wind of this duration has had time to exploit.

    Inverts the Carter (1982) relation  t~ = 68.8 x~^0.67.
    """
    t_tilde = dimensionless_duration(duration, U10)
    if t_tilde <= 0:
        return 0.0
    x_tilde_eff = (t_tilde / 68.8) ** (1.0 / 0.67)
    return x_tilde_eff * U10**2 / G


def sea_state(U10, fetch, duration=np.inf):
    """Work out which of the three growth limits is binding, and the resulting
    peak frequency and scale constant.

    A sea can be limited by three things: how far the wind blows over water
    (fetch), how long it blows (duration), or nothing at all, in which case it
    saturates at the Pierson-Moskowitz state and stops growing no matter what.
    Whichever gives the smallest effective fetch wins.
    """
    if U10 <= 0:
        raise ValueError("U10 must be positive")

    x_fetch = dimensionless_fetch(fetch, U10)
    x_dur = (dimensionless_fetch(duration_limited_fetch(duration, U10), U10)
             if np.isfinite(duration) else np.inf)

    x_tilde = min(x_fetch, x_dur, FULLY_DEVELOPED_FETCH)

    if x_tilde == FULLY_DEVELOPED_FETCH:
        regime = "fully-developed"
    elif x_dur < x_fetch:
        regime = "duration-limited"
    else:
        regime = "fetch-limited"

    f_p = 3.5 * (G / U10) * x_tilde ** (-0.33)
    alpha = 0.076 * x_tilde ** (-0.22)

    return SeaState(U10=U10, fetch=fetch, duration=duration, x_tilde=x_tilde,
                    f_p=f_p, alpha=alpha, regime=regime)


def significant_height_scaling(U10, fetch, duration=np.inf):
    """H_s straight from the fetch law, bypassing the spectrum entirely.

        g H_s / U10^2 = 0.0016 sqrt(x~)

    Useful precisely because it is independent: if integrating the JONSWAP
    spectrum does not reproduce this to within a percent or two, something is
    wrong with the integration. See tests/test_spectra.py.
    """
    s = sea_state(U10, fetch, duration)
    return 0.0016 * np.sqrt(s.x_tilde) * U10**2 / G
