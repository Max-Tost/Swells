"""Turning a spectrum back into a sea surface.

Derived in lessons/07-why-waves-come-in-sets.md.

A spectrum tells you how much variance sits at each frequency, and nothing at
all about phase. To watch sets arrive you need an actual time series, so we put
the phases back in at random. The result is one realisation of a sea with that
spectrum -- not *the* sea, but a perfectly good sample of it.
"""

import numpy as np

from .util import trapz


def surface_elevation(f, S, duration, dt, seed=None, return_components=False):
    """Synthesise eta(t) from a variance density spectrum by random phases.

    eta(t) = sum_j a_j cos(2 pi f_j t + phi_j),   a_j = sqrt(2 S(f_j) df)

    with phi_j drawn uniformly. The a_j come out right because the variance of
    a single cosine of amplitude a is a^2/2, and we want each band df to
    contribute S df.

    By the central limit theorem the sum of many such components is Gaussian,
    which is the assumption every statistic in lesson 07 rests on.

    Built with an inverse FFT, so it is O(N log N) rather than O(N * n_freq).

    Returns (t, eta), or (t, eta, dict) with the amplitudes and phases.
    """
    n = int(round(duration / dt))
    if n < 4:
        raise ValueError("duration too short for this dt")

    df = 1.0 / (n * dt)
    f_grid = np.fft.rfftfreq(n, dt)

    S_grid = np.interp(f_grid, np.asarray(f), np.asarray(S), left=0.0, right=0.0)
    S_grid[0] = 0.0                      # no mean offset; the sea is at z = 0
    a = np.sqrt(2.0 * np.clip(S_grid, 0.0, None) * df)

    rng = np.random.default_rng(seed)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=f_grid.size)
    phi[0] = 0.0

    # irfft normalises by n, so scale by n/2 to make bin j come out with
    # amplitude a_j.
    X = (n / 2.0) * a * np.exp(1j * phi)
    if n % 2 == 0:
        X[-1] = X[-1].real     # Nyquist bin must be real
    eta = np.fft.irfft(X, n=n)

    t = np.arange(n) * dt
    if return_components:
        return t, eta, {"f": f_grid, "a": a, "phi": phi}
    return t, eta


def synthesise_from_spectrogram(f, S_ft, t_spec, duration, dt, seed=None):
    """Synthesise a record from a slowly evolving spectrum.

    Blends the spectrum across the record so a multi-day swell event can be
    played out with the period falling as the wind sea catches up. Phases stay
    fixed across the blend so the record is continuous.
    """
    n = int(round(duration / dt))
    df = 1.0 / (n * dt)
    f_grid = np.fft.rfftfreq(n, dt)
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=f_grid.size)

    t = np.arange(n) * dt

    # Regrid the spectrogram onto (t, f_grid): first in frequency, then in time.
    S_on_fgrid = np.array([np.interp(f_grid, f, S_ft[i]) for i in range(len(t_spec))])
    S_interp = np.array([np.interp(t, t_spec, S_on_fgrid[:, j])
                         for j in range(f_grid.size)]).T
    S_interp = np.clip(S_interp, 0.0, None)

    a_of_t = np.sqrt(2.0 * S_interp * df)
    eta = np.sum(a_of_t * np.cos(2.0 * np.pi * f_grid[None, :] * t[:, None]
                                 + phi[None, :]), axis=1)
    return t, eta


def check_variance(f, S, eta):
    """m_0 from the spectrum versus the variance of the record.

    These agree to within sampling error for a long enough record. If they do
    not, the synthesis is wrong. See tests/test_synth.py.
    """
    m0_spec = float(trapz(S, f))
    return m0_spec, float(np.var(eta))
