"""End to end: the whole chain must hang together and behave sensibly."""

import numpy as np
import pytest

from swells import Coast, Storm, run, surf_report
from swells.nonlinear import (BF_CRITICAL_KH, benjamin_feir_growth_rate,
                              benjamin_feir_index, is_focusing,
                              nls_coefficients, unstable_sideband_range)


def test_default_run_is_physical():
    r = run(record_hours=0.5)
    assert 0.5 < r.peak["Hm0"] < 6.0
    assert 8.0 < r.peak["Tp"] < 25.0
    assert r.breaking is not None
    assert 1.0 < r.grouping["mean_run_length"] < 10.0
    assert 30.0 < r.grouping["set_interval"] < 1200.0


def test_report_renders():
    text = surf_report(run(record_hours=0.25))
    for section in ["GENERATION", "CROSSING", "AT THE BUOY", "SETS", "SURF"]:
        assert section in text


def test_distance_makes_swell_smaller_cleaner_and_longer_lasting():
    """The three signatures of a long journey."""
    near = run(storm=Storm(distance_km=1000), record_hours=0.25)
    far = run(storm=Storm(distance_km=12000), record_hours=0.25)
    assert far.peak["Hm0"] < near.peak["Hm0"]
    assert far.peak["nu"] < near.peak["nu"]
    assert far.grouping["mean_run_length"] > near.grouping["mean_run_length"]


def test_stronger_wind_makes_bigger_longer_waves():
    weak = run(storm=Storm(U10=12.0), record_hours=0.25)
    strong = run(storm=Storm(U10=30.0), record_hours=0.25)
    assert strong.peak["Hm0"] > weak.peak["Hm0"]
    assert strong.peak["Tp"] > weak.peak["Tp"]


def test_longer_fetch_makes_longer_period():
    short = run(storm=Storm(fetch_km=100), record_hours=0.25)
    long = run(storm=Storm(fetch_km=1500), record_hours=0.25)
    assert long.peak["Tp"] > short.peak["Tp"]


def test_peak_arrives_when_the_dispersion_law_says_it_should():
    r = run(storm=Storm(distance_km=6000), record_hours=0.25)
    from swells.propagate import travel_time
    expected_h = travel_time(1 / r.peak["Tp"], r.storm.distance) / 3600.0
    assert r.peak["time_h"] == pytest.approx(expected_h, rel=0.10)


def test_steep_beach_changes_breaker_type():
    mellow = run(coast=Coast(slope=0.01), record_hours=0.25)
    steep = run(coast=Coast(slope=0.10), record_hours=0.25)
    assert steep.breaking.xi0 > mellow.breaking.xi0


# --- nonlinear module ------------------------------------------------------

def test_nls_deep_water_coefficients():
    omega = 2 * np.pi / 12.0
    alpha, beta, k = nls_coefficients(omega)
    assert alpha == pytest.approx(-omega / (8 * k**2))
    assert beta == pytest.approx(-0.5 * omega * k**2)


def test_focusing_switches_off_in_shallow_water():
    k = 0.05
    assert is_focusing(k, 100.0)
    assert not is_focusing(k, BF_CRITICAL_KH / k * 0.9)


def test_benjamin_feir_growth_timescale():
    """Gamma = omega eps^2 / 2, so the e-folding time is 1/Gamma = 2/(omega eps^2)
    seconds, which is 1/(pi eps^2) wave periods. At eps = 0.1 that is 31.8
    periods -- about eight minutes for a 15 s wave. At a realistic swell
    steepness of 0.02 it is 800 periods, several hours."""
    omega = 2 * np.pi / 15.0
    k = omega**2 / 9.80665
    a = 0.1 / k
    gamma = benjamin_feir_growth_rate(omega, k, a)
    n_periods = (1.0 / gamma) / (2 * np.pi / omega)
    assert n_periods == pytest.approx(1 / (np.pi * 0.1**2), rel=1e-9)
    assert n_periods == pytest.approx(31.83, abs=0.01)


def test_unstable_band():
    k, a = 0.028, 1.0
    K_max, K_cut = unstable_sideband_range(k, a)
    assert K_cut == pytest.approx(np.sqrt(2) * K_max)


def test_bfi_rises_as_the_spectrum_narrows():
    assert benjamin_feir_index(0.03, 1.0, 0.05) > \
        benjamin_feir_index(0.03, 1.0, 0.30)
