"""Transoceanic propagation: the dispersive sorting of a storm spectrum.

Derived in lessons/05-the-great-sorting.md and lessons/06-what-survives.md.

Once swell leaves the storm there is no wind input and essentially no
dissipation, so the physics reduces to bookkeeping: which frequency arrives
when, and how much of it survives the journey.
"""

import numpy as np

from .constants import G, R_EARTH
from .dispersion import deep_group_speed
from .util import trapz


# ---------------------------------------------------------------------------
# The kinematics of arrival
# ---------------------------------------------------------------------------

def travel_time(f, D):
    """t = D / c_g(f) = 4 pi f D / g.

    Linear in f, which is the entire reason a spectrogram of arriving swell
    shows a straight ridge.
    """
    return 4.0 * np.pi * np.asarray(f, dtype=np.float64) * D / G


def arrival_frequency(t, D):
    """The inverse: f(t) = g t / (4 pi D). Frequency rises with time."""
    return G * np.asarray(t, dtype=np.float64) / (4.0 * np.pi * D)


def chirp_rate(D):
    """df/dt = g / (4 pi D), in Hz/s. The slope of the ridge."""
    return G / (4.0 * np.pi * D)


def distance_from_chirp(df_dt):
    """The inverse problem: read the storm's distance off a buoy spectrogram.

    D = g / (4 pi df/dt).
    """
    return G / (4.0 * np.pi * df_dt)


def origin_time_from_ridge(df_dt, f_intercept_time):
    """Given the ridge slope and the time at which it extrapolates to f = 0,
    that intercept *is* the moment the storm radiated. No extra work needed.
    """
    return f_intercept_time


# ---------------------------------------------------------------------------
# What survives
# ---------------------------------------------------------------------------

def geometric_spreading(D, D_ref, R=R_EARTH):
    """Energy attenuation factor from spreading along great circles.

    A fan of rays leaving the storm with azimuthal width dtheta occupies an arc
    of length R sin(D/R) dtheta at range D. Energy flux through that arc is
    conserved, so E ~ 1 / (R sin(D/R)).

    Returns E(D)/E(D_ref), so amplitudes scale as its square root.

    For D << R this is the flat-Earth D_ref/D. Past a quarter of the way round
    the planet, sin(D/R) starts *decreasing* and the swell re-focuses. That is
    real: it is why the antipode of a Southern Ocean storm can get more swell
    than somewhere closer.
    """
    D = np.asarray(D, dtype=np.float64)
    arc = R * np.sin(np.clip(D / R, 1e-9, np.pi - 1e-9))
    arc_ref = R * np.sin(np.clip(D_ref / R, 1e-9, np.pi - 1e-9))
    return arc_ref / arc


def arrival_bandwidth(f, D, storm_duration, fetch):
    """Instantaneous frequency bandwidth arriving at range D.

        df ~ (g/4 pi D) * T_storm  +  f * (fetch / D)

    First term: the storm radiated for a finite time, so each frequency arrives
    smeared over that time, which the f-t mapping converts into a frequency
    width. Second term: the storm has finite size, so its near and far edges
    are at different ranges.

    The storm's extent is taken to be the fetch itself -- the wind blows across
    the storm, so the patch radiating at you is one fetch wide. That ties the
    only two lengths in the source problem together instead of leaving them
    free to contradict each other; see lesson 05 for what it assumes.

    Both terms shrink as 1/D. This is why distant swell is clean.
    """
    f = np.asarray(f, dtype=np.float64)
    return chirp_rate(D) * storm_duration + f * (fetch / D)


def arrival_duration(f, storm_duration, fetch):
    """How long a single frequency keeps arriving, in seconds.

    Dispersion does not stretch a single frequency component -- they all travel
    at the same c_g. The stretching comes only from the source: the storm blew
    for T_storm, and it is one fetch across.
    """
    return storm_duration + fetch / deep_group_speed(f)


def spectrogram(f, S_source, D, storm_duration, fetch, t,
                R=R_EARTH, D_ref=None):
    """The spectrum a buoy at range D sees, as a function of time.

    S(f, t) = S_source(f) * spreading(D) * (T_storm / dT_eff(f)) * w(t - t_a(f))

    The three factors after the source spectrum are, in order: energy lost to
    geometric spreading; the drop in spectral *density* because the same
    variance now arrives spread over a longer window than it was radiated over;
    and a normalised window placing it at the right arrival time.

    This is a kinematic model. It transports energy correctly and puts it in
    the right place at the right time, but it does not solve the action balance
    equation -- there is no nonlinear transfer and no dissipation in transit.
    Over deep ocean those really are small (lesson 06), which is what makes the
    shortcut defensible.

    Returns an array of shape (len(t), len(f)).
    """
    f = np.asarray(f, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    S_source = np.asarray(S_source, dtype=np.float64)

    if D_ref is None:
        D_ref = max(0.5 * fetch, 1.0)   # centre of the storm to its edge

    spread = geometric_spreading(D, D_ref, R)

    t_a = travel_time(f, D)                                  # (nf,)
    dT = arrival_duration(f, storm_duration, fetch)          # (nf,)
    stretch = storm_duration / dT

    # Gaussian window whose full width at ~2 sigma is dT.
    sigma_t = np.maximum(dT / 4.0, 1e-6)
    w = np.exp(-0.5 * ((t[:, None] - t_a[None, :]) / sigma_t[None, :]) ** 2)

    return S_source[None, :] * spread * stretch[None, :] * w


def spectrogram_ridge(f, S_ft):
    """The peak frequency at each time -- the ridge you would fit by eye."""
    return np.asarray(f)[np.argmax(S_ft, axis=1)]


def fit_chirp(t, f_ridge, weights=None):
    """Least-squares fit f = a t + b to a ridge. Returns (df_dt, intercept_time).

    Recovering `a` gives the storm distance and the zero crossing gives the
    time it blew. This is the inverse problem worked in lesson 05.
    """
    t = np.asarray(t, dtype=np.float64)
    f_ridge = np.asarray(f_ridge, dtype=np.float64)
    w = np.ones_like(t) if weights is None else np.asarray(weights, dtype=np.float64)
    ok = w > 0
    a, b = np.polyfit(t[ok], f_ridge[ok], 1, w=np.sqrt(w[ok]))
    return float(a), float(-b / a) if a != 0 else np.nan


def sea_state_at(f, S_ft):
    """H_m0 and peak period as functions of time, from a spectrogram."""
    m0 = trapz(S_ft, f, axis=1)
    Hm0 = 4.0 * np.sqrt(np.clip(m0, 0.0, None))
    f_pk = np.asarray(f)[np.argmax(S_ft, axis=1)]
    with np.errstate(divide="ignore"):
        Tp = np.where(f_pk > 0, 1.0 / f_pk, 0.0)
    return Hm0, Tp
