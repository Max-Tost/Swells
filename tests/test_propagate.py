"""The ocean crossing: arrival times, the chirp, and the inverse problem."""

import numpy as np
import pytest

from swells.propagate import (arrival_bandwidth, arrival_frequency, chirp_rate,
                              distance_from_chirp, fit_chirp,
                              geometric_spreading, sea_state_at, spectrogram,
                              spectrogram_ridge, travel_time)
from swells.spectra import Hm0, storm_spectrum

F = np.linspace(0.02, 0.25, 900)


def test_benchmark_2_transit():
    """D = 4000 km, T_p = 14.918 s."""
    f_p = 0.067034
    D = 4.0e6
    assert travel_time(f_p, D) / 3600.0 == pytest.approx(95.44, abs=0.02)
    assert chirp_rate(D) == pytest.approx(1.95097e-7, rel=1e-4)
    assert chirp_rate(D) * 86400 == pytest.approx(0.016856, abs=1e-5)


def test_arrival_frequency_inverts_travel_time():
    D = 7.5e6
    f = np.linspace(0.03, 0.2, 50)
    assert np.allclose(arrival_frequency(travel_time(f, D), D), f)


def test_distance_from_chirp_inverts_chirp_rate():
    for D in [1e6, 4e6, 1.2e7]:
        assert distance_from_chirp(chirp_rate(D)) == pytest.approx(D, rel=1e-12)


def test_the_inverse_problem_recovers_the_storm():
    """Synthesise a spectrogram for a known storm, then fit the ridge and see
    whether we get the distance back. This is the punchline of lesson 05."""
    D_true = 8.0e6
    S0, sea = storm_spectrum(F, 22.0, 800e3, 36 * 3600.0)
    t = np.linspace(0.4 * travel_time(F.min(), D_true),
                    1.2 * travel_time(F.max(), D_true), 1200)
    S_ft = spectrogram(F, S0, D_true, 36 * 3600.0, 800e3, t)

    Hm0_t, _ = sea_state_at(F, S_ft)
    ridge = spectrogram_ridge(F, S_ft)
    df_dt, t0 = fit_chirp(t, ridge, weights=Hm0_t**2)

    assert distance_from_chirp(df_dt) == pytest.approx(D_true, rel=0.05)
    assert abs(t0) < 0.05 * travel_time(sea.f_p, D_true)


def test_swell_decays_with_distance():
    """The bug in the draft: arriving Hm0 was independent of D."""
    S0, _ = storm_spectrum(F, 25.0, 600e3, 36 * 3600.0)
    peaks = []
    for D in [1e6, 4e6, 1.2e7]:
        t = np.linspace(0.5 * travel_time(F.min(), D),
                        1.3 * travel_time(F.max(), D), 800)
        S_ft = spectrogram(F, S0, D, 36 * 3600.0, 600e3, t)
        peaks.append(Hm0(F, S_ft[np.argmax(S_ft.max(axis=1))]))
    assert peaks[0] > peaks[1] > peaks[2]


def test_geometric_spreading_is_nearly_inverse_distance_when_flat():
    """Close in, E ~ 1/D. Not exactly: the arc R sin(D/R) is shorter than the
    flat-Earth D, so a sphere spreads energy slightly less than a plane. At
    2000 km that is a 1.3% effect, and it grows from there."""
    assert geometric_spreading(1e6, 1e6) == pytest.approx(1.0)
    assert geometric_spreading(2e6, 1e6) == pytest.approx(0.5, rel=0.02)
    assert geometric_spreading(2e6, 1e6) > 0.5


def test_antipodal_refocusing():
    """Past a quarter of the way round the Earth the arc shrinks again."""
    R = 6.371e6
    quarter = 0.5 * np.pi * R
    assert geometric_spreading(quarter + 2e6, 1e6) > \
        geometric_spreading(quarter, 1e6)


def test_bandwidth_narrows_with_distance():
    """Why distant swell is clean."""
    f = 0.07
    wide = arrival_bandwidth(f, 1e6, 36 * 3600.0, 800e3)
    narrow = arrival_bandwidth(f, 1.2e7, 36 * 3600.0, 800e3)
    assert narrow < wide / 5


def test_long_period_arrives_first():
    D = 6e6
    assert travel_time(1 / 20.0, D) < travel_time(1 / 10.0, D)
    # and exactly twice as fast, since c_g ~ 1/f
    assert travel_time(1 / 10.0, D) == pytest.approx(2 * travel_time(1 / 20.0, D))
