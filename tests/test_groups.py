"""Sets: the theory must agree with what you get by counting."""

import numpy as np

from swells.util import trapz
import pytest

from swells.groups import (beat_period, bivariate_rayleigh_pdf, count_runs,
                           envelope, envelope_correlation, group_statistics,
                           height_statistics, mean_run_length,
                           measured_group_statistics, p22,
                           rayleigh_exceedance, zero_crossing_waves)
from swells.spectra import Hm0, bandwidth_nu, storm_spectrum
from swells.synth import check_variance, surface_elevation

F = np.linspace(0.02, 0.40, 1500)


def gaussian_band(f0, width):
    S = np.exp(-((F - f0) ** 2) / (2 * width**2))
    return S / trapz(S, F) * 1.0        # m0 = 1 m^2, so Hm0 = 4 m


# --- the bivariate model ---------------------------------------------------

def test_bivariate_rayleigh_normalises():
    x = np.linspace(1e-6, 8.0, 900)
    X1, X2 = np.meshgrid(x, x, indexing="ij")
    for kappa in [0.0, 0.5, 0.9]:
        p = bivariate_rayleigh_pdf(X1, X2, kappa)
        total = trapz(trapz(p, x, axis=1), x)
        assert total == pytest.approx(1.0, rel=2e-3)


def test_bivariate_factorises_at_zero_correlation():
    x = np.linspace(1e-6, 6.0, 300)
    X1, X2 = np.meshgrid(x, x, indexing="ij")
    joint = bivariate_rayleigh_pdf(X1, X2, 0.0)
    marg = 2 * x * np.exp(-(x**2))
    assert np.allclose(joint, np.outer(marg, marg), atol=1e-10)


def test_p22_at_zero_correlation_is_the_unconditional_probability():
    """Independent waves: knowing this one was big tells you nothing.
    p22 = exp(-xc^2) = e^-2 = 0.135, giving 1.16 waves per 'set' -- i.e. no
    grouping at all. This is the null hypothesis every set has to beat."""
    assert p22(0.0) == pytest.approx(np.exp(-2.0), rel=2e-3)
    assert mean_run_length(p22(0.0)) == pytest.approx(1.156, rel=5e-3)


def test_p22_increases_with_correlation():
    ps = [p22(k) for k in [0.0, 0.3, 0.6, 0.8, 0.95]]
    assert all(b > a for a, b in zip(ps, ps[1:]))


def test_narrow_spectra_are_groupier():
    """The whole claim of lesson 07 in one assertion."""
    narrow = group_statistics(F, gaussian_band(0.07, 0.003))
    broad = group_statistics(F, gaussian_band(0.07, 0.02))
    assert narrow["kappa"] > broad["kappa"]
    assert narrow["mean_run_length"] > broad["mean_run_length"]


def test_swell_is_groupier_than_wind_sea():
    swell = group_statistics(F, gaussian_band(0.06, 0.004))
    sea, _ = storm_spectrum(F, 25.0, 300e3)
    windsea = group_statistics(F, sea)
    assert swell["mean_run_length"] > windsea["mean_run_length"]
    assert bandwidth_nu(F, gaussian_band(0.06, 0.004)) < bandwidth_nu(F, sea)


# --- Rayleigh height statistics -------------------------------------------

def test_rayleigh_ratios():
    st = height_statistics(4.0)
    assert st["H_rms"] / st["Hm0"] == pytest.approx(1 / np.sqrt(2), rel=1e-3)
    assert st["H_mean"] / st["Hm0"] == pytest.approx(0.626, abs=2e-3)
    assert st["H_tenth"] / st["Hm0"] == pytest.approx(1.272, abs=2e-3)
    assert st["H_hundredth"] / st["Hm0"] == pytest.approx(1.668, abs=2e-3)


def test_rayleigh_exceedance_at_hs():
    """One wave in eight exceeds H_s."""
    assert rayleigh_exceedance(4.0, 4.0) == pytest.approx(np.exp(-2), rel=1e-9)


# --- synthesis -------------------------------------------------------------

def test_synthesis_reproduces_m0():
    S = gaussian_band(0.07, 0.01)
    t, eta = surface_elevation(F, S, 6 * 3600.0, 0.5, seed=1)
    m0_spec, var = check_variance(F, S, eta)
    assert var == pytest.approx(m0_spec, rel=0.02)
    assert 4 * np.std(eta) == pytest.approx(Hm0(F, S), rel=0.02)


def test_synthesised_record_is_gaussian():
    S = gaussian_band(0.07, 0.015)
    _, eta = surface_elevation(F, S, 12 * 3600.0, 0.5, seed=2)
    from scipy.stats import kurtosis, skew
    assert abs(skew(eta)) < 0.05
    assert abs(kurtosis(eta)) < 0.1


def test_zero_crossing_finds_the_right_number_of_waves():
    """A pure 10 s sine over an hour must give 360 waves of height 2 m."""
    t = np.arange(0, 3600, 0.1)
    eta = 1.0 * np.cos(2 * np.pi * t / 10.0)
    H, T = zero_crossing_waves(t, eta)
    assert H.size == pytest.approx(359, abs=2)
    assert np.median(H) == pytest.approx(2.0, rel=1e-3)
    assert np.median(T) == pytest.approx(10.0, rel=1e-3)


def test_envelope_of_a_beat_is_the_beat():
    """Two equal sines at f1, f2 give an envelope 2a|cos(pi(f1-f2)t)|."""
    f1, f2 = 1 / 14.0, 1 / 15.0
    t = np.arange(0, 4000, 0.05)
    eta = np.cos(2 * np.pi * f1 * t) + np.cos(2 * np.pi * f2 * t)
    env = envelope(eta)
    inner = (t > 500) & (t < 3500)
    assert env[inner].max() == pytest.approx(2.0, rel=0.02)
    assert env[inner].min() == pytest.approx(0.0, abs=0.05)
    assert beat_period(f1, f2) == pytest.approx(210.0, rel=1e-6)


# --- the two routes must agree --------------------------------------------

def test_predicted_and_counted_run_lengths_agree():
    """Theory from the spectrum vs counting in a synthesised record.

    This is the test that matters. The draft's p22 was an ad-hoc fudge; this
    catches that.
    """
    for width in [0.004, 0.010, 0.020]:
        S = gaussian_band(0.07, width)
        predicted = group_statistics(F, S)["mean_run_length"]
        t, eta = surface_elevation(F, S, 40 * 3600.0, 0.5, seed=7)
        measured = measured_group_statistics(t, eta)["mean_run_length"]
        assert measured == pytest.approx(predicted, rel=0.15), (
            f"width={width}: predicted {predicted:.2f}, counted {measured:.2f}")


def test_count_runs():
    h = np.array([1, 3, 3, 1, 1, 3, 1, 3, 3, 3])
    assert list(count_runs(h, 2)) == [2, 1, 3]
