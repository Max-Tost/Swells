"""Wave spectra: JONSWAP, Pierson-Moskowitz, moments, directional spreading.

Derived in lessons/04-the-spectrum.md.

S(f) is a variance density, in m^2/Hz: integrate it over a band of frequencies
and you get the contribution of that band to the variance of the sea surface.
Every wave-height statistic in this project is a moment of S.
"""

from dataclasses import replace

import numpy as np
from scipy.special import gammaln

from .constants import G
from .fetch import sea_state
from .util import trapz


def pierson_moskowitz(f, U10=None, f_p=None, alpha=0.0081):
    """The fully developed sea (Pierson & Moskowitz 1964).

    Either give U10 and let f_p = 0.13 g / U10, or give f_p directly.
    """
    f = np.asarray(f, dtype=np.float64)
    if f_p is None:
        if U10 is None:
            raise ValueError("give either U10 or f_p")
        f_p = 0.13 * G / U10
    return _base_spectrum(f, alpha, f_p)


def jonswap(f, alpha, f_p, gamma=3.3, sigma_a=0.07, sigma_b=0.09):
    """The fetch-limited spectrum (Hasselmann et al. 1973).

    Pierson-Moskowitz multiplied by a Gaussian bump of height gamma sitting on
    the peak. A young, actively driven sea has a sharper peak than a mature
    one: gamma is the sharpness.

    gamma = 3.3 is the JONSWAP mean. Individual storms in that dataset ranged
    from about 1 to 7, so treat 3.3 as a central value, not a constant of
    nature. gamma = 1 recovers Pierson-Moskowitz exactly.
    """
    f = np.asarray(f, dtype=np.float64)
    sigma = np.where(f <= f_p, sigma_a, sigma_b)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.exp(-((f - f_p) ** 2) / (2.0 * sigma**2 * f_p**2))
    return _base_spectrum(f, alpha, f_p) * gamma**r


def _base_spectrum(f, alpha, f_p):
    """alpha g^2 (2 pi)^-4 f^-5 exp[-1.25 (f_p/f)^4], with f = 0 handled."""
    f = np.asarray(f, dtype=np.float64)
    out = np.zeros_like(f)
    ok = f > 0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        out[ok] = (alpha * G**2 / ((2.0 * np.pi) ** 4 * f[ok] ** 5)
                   * np.exp(-1.25 * (f_p / f[ok]) ** 4))
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def storm_spectrum(f, U10, fetch, duration=np.inf, gamma=3.3, calibrate=True):
    """The spectrum a storm of given wind, fetch and duration produces.

    With calibrate=True (the default) the JONSWAP shape is rescaled so that its
    total energy obeys the JONSWAP energy law rather than whatever the alpha
    and f_p fits happen to imply. Those two disagree by up to 25%; see the long
    comment in fetch.py for why, and why we side with the energy law.

    The scale factor is computed on an internal wide frequency grid, so it does
    not depend on the grid you pass in. Returns (S, SeaState); the factor used
    is on the SeaState as .energy_scale.
    """
    s = sea_state(U10, fetch, duration)
    scale = 1.0
    if calibrate:
        f_ref = np.logspace(np.log10(0.15 * s.f_p), np.log10(25.0 * s.f_p), 4000)
        m0_raw = moment(f_ref, jonswap(f_ref, s.alpha, s.f_p, gamma=gamma), 0)
        if m0_raw > 0:
            scale = s.target_m0 / m0_raw
    s = replace(s, energy_scale=scale)
    return jonswap(f, s.alpha, s.f_p, gamma=gamma) * scale, s


# ---------------------------------------------------------------------------
# Moments and the statistics built from them
# ---------------------------------------------------------------------------

def moment(f, S, n=0):
    """m_n = int f^n S(f) df."""
    f = np.asarray(f, dtype=np.float64)
    return float(trapz(np.asarray(S, dtype=np.float64) * f**n, f))


def moments(f, S, orders=(0, 1, 2, 4)):
    return {n: moment(f, S, n) for n in orders}


def Hm0(f, S):
    """H_m0 = 4 sqrt(m_0). The spectral significant wave height."""
    return 4.0 * np.sqrt(max(moment(f, S, 0), 0.0))


def peak_frequency(f, S):
    return float(np.asarray(f)[int(np.argmax(S))])


def mean_frequency(f, S):
    """f_bar = m_1/m_0, the spectral centroid. Its reciprocal is the mean period."""
    m0 = moment(f, S, 0)
    return moment(f, S, 1) / m0 if m0 > 0 else 0.0


def mean_period(f, S):
    """T_m01 = m_0/m_1. This is the period that counts waves."""
    fb = mean_frequency(f, S)
    return 1.0 / fb if fb > 0 else 0.0


def bandwidth_nu(f, S, f_cut=None):
    """The narrowness parameter nu = sqrt(m0 m2 / m1^2 - 1).

    Zero for a pure sine wave, order 0.3 for a wind sea, order 0.03 for clean
    distant swell. This single number is what decides how groupy the sea is,
    which is the whole content of lessons/07.

    A warning that matters. m_2 weights the integrand by f^2, and a JONSWAP
    tail goes as f^-5, so the m_2 integrand decays only as f^-3. The integral
    converges, but slowly, and nu is therefore sensitive to where you truncate.
    Worse, since a JONSWAP is nearly self-similar in f/f_p, what nu really
    measures is the ratio f_max/f_p -- so an old sea with a low f_p reports a
    *larger* nu than a young one on the same grid, which is the opposite of the
    physical intent.

    Pass f_cut to truncate at a fixed multiple of the peak (3 f_p is the usual
    choice) and get a number that means what you wanted it to mean.
    """
    f = np.asarray(f, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)
    if f_cut is not None:
        keep = f <= f_cut
        f, S = f[keep], S[keep]
    m0, m1, m2 = moment(f, S, 0), moment(f, S, 1), moment(f, S, 2)
    if m0 <= 0 or m1 <= 0:
        return 0.0
    return float(np.sqrt(max(m0 * m2 / m1**2 - 1.0, 0.0)))


# ---------------------------------------------------------------------------
# Directional spreading
# ---------------------------------------------------------------------------

def spreading_exponent(f, f_p, s_max=10.0):
    """Mitsuyasu et al. (1975): s peaks at f_p and falls off either side.

    s_max ~ 10 for a wind sea under active forcing, 25-75 for old swell that
    has had the whole ocean to sort itself out directionally.
    """
    f = np.asarray(f, dtype=np.float64)
    r = np.where(f <= f_p, 5.0, -2.5)
    return s_max * (f / f_p) ** r


def cos2s_spreading(theta, theta_mean, s):
    """D(theta) proportional to cos^{2s}((theta - theta_mean)/2), normalised so
    that the integral over all directions is 1.

    The normalisation is 2^(2s-1) Gamma(s+1)^2 / (pi Gamma(2s+1)), computed
    through log-gamma because s can reach 75 and the factorials overflow.
    """
    theta = np.asarray(theta, dtype=np.float64)
    s = np.asarray(s, dtype=np.float64)
    log_norm = ((2.0 * s - 1.0) * np.log(2.0) + 2.0 * gammaln(s + 1.0)
                - np.log(np.pi) - gammaln(2.0 * s + 1.0))
    half = 0.5 * (theta - theta_mean)
    c = np.cos(half)
    return np.where(c > 0, np.exp(log_norm) * np.abs(c) ** (2.0 * s), 0.0)


def directional_spectrum(f, theta, S, f_p, theta_mean=0.0, s_max=10.0):
    """S(f, theta) = S(f) D(f, theta), shaped (len(f), len(theta))."""
    s = spreading_exponent(f, f_p, s_max)
    D = cos2s_spreading(theta[None, :], theta_mean, s[:, None])
    return np.asarray(S, dtype=np.float64)[:, None] * D
