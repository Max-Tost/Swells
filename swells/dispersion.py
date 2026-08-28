"""Linear surface gravity waves: the dispersion relation and what follows from it.

Derived in lessons/01-the-water-wave-problem.md (the relation itself) and
lessons/02-group-velocity.md (phase speed, group speed, energy, energy flux).

Everything here is linear theory: amplitudes small compared with wavelength,
inviscid, irrotational, flat bottom locally.
"""

import numpy as np

from .constants import G, RHO_W


def wavenumber(omega, h, tol=1e-13, max_iter=100):
    """Solve  omega^2 = g k tanh(k h)  for the wavenumber k.

    There is no closed form, so we use Newton-Raphson from Eckart's
    approximation, which is good to a few percent everywhere and converges in
    three or four iterations.

    omega : angular frequency, rad/s
    h     : water depth, m
    Broadcasts over both arguments.
    """
    omega = np.asarray(omega, dtype=np.float64)
    h = np.asarray(h, dtype=np.float64)
    omega, h = np.broadcast_arrays(omega, h)

    k_deep = omega**2 / G                       # the deep-water answer
    # Eckart (1952): k ~ k_deep / sqrt(tanh(k_deep h)). Exact in both limits.
    k = k_deep / np.sqrt(np.tanh(np.clip(k_deep * h, 1e-12, 700.0)))

    for _ in range(max_iter):
        kh = np.clip(k * h, 1e-12, 700.0)
        t = np.tanh(kh)
        f = G * k * t - omega**2
        df = G * t + G * k * h * (1.0 - t**2)
        step = f / df
        k = k - step
        if np.all(np.abs(step) <= tol * np.abs(k)):
            break

    return k


def wavelength(omega, h):
    return 2.0 * np.pi / wavenumber(omega, h)


def phase_speed(omega, h):
    """c = omega / k. The speed of an individual crest."""
    return omega / wavenumber(omega, h)


def n_ratio(omega, h):
    """n = c_g / c = (1/2)(1 + 2kh/sinh 2kh).

    Runs from 1/2 in deep water to 1 in shallow water.
    """
    kh = np.clip(wavenumber(omega, h) * h, 1e-12, 350.0)
    return 0.5 * (1.0 + 2.0 * kh / np.sinh(2.0 * kh))


def group_speed(omega, h):
    """c_g = dω/dk = n c. The speed at which energy travels."""
    return n_ratio(omega, h) * phase_speed(omega, h)


# ---------------------------------------------------------------------------
# Deep-water shortcuts, in terms of ordinary frequency f = 1/T.
# These are the workhorses for the open-ocean part of the problem, where
# kh >> 1 and the tanh is 1 to within a rounding error.
# ---------------------------------------------------------------------------

def deep_phase_speed(f):
    """c = g / (2 pi f) = g T / (2 pi).  A 20 s wave: 31.2 m/s."""
    return G / (2.0 * np.pi * np.asarray(f, dtype=np.float64))


def deep_group_speed(f):
    """c_g = g / (4 pi f) = g T / (4 pi).  Half the phase speed.

    Sanity anchor: T = 20 s gives c_g = 9.80665*20/(4*pi) = 15.61 m/s
    = 56.2 km/h, so 10,000 km takes 7.4 days.
    """
    return G / (4.0 * np.pi * np.asarray(f, dtype=np.float64))


def deep_wavelength(f):
    """L = g T^2 / (2 pi).  A 15 s wave is 351 m long."""
    f = np.asarray(f, dtype=np.float64)
    return G / (2.0 * np.pi * f**2)


def is_deep(f, h):
    """Deep water means kh > pi, i.e. h > L/2. Returns a boolean array."""
    omega = 2.0 * np.pi * np.asarray(f, dtype=np.float64)
    return wavenumber(omega, h) * h > np.pi


def is_shallow(f, h):
    """Shallow water means kh < pi/10, i.e. h < L/20."""
    omega = 2.0 * np.pi * np.asarray(f, dtype=np.float64)
    return wavenumber(omega, h) * h < np.pi / 10.0


# ---------------------------------------------------------------------------
# Energy
# ---------------------------------------------------------------------------

def energy_density(H, rho=RHO_W):
    """E = (1/8) rho g H^2, in J/m^2. Half potential, half kinetic."""
    return 0.125 * rho * G * np.asarray(H, dtype=np.float64) ** 2


def energy_flux(H, omega, h, rho=RHO_W):
    """F = E c_g, in W/m of crest. This is the quantity conserved in shoaling."""
    return energy_density(H, rho) * group_speed(omega, h)
