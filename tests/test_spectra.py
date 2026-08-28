"""The spectrum must reproduce the fetch laws it was built from."""

import numpy as np

from swells.util import trapz
import pytest

from swells.constants import G
from swells.fetch import (FULLY_DEVELOPED_FETCH, dimensionless_fetch,
                          duration_limited_fetch, sea_state,
                          significant_height_scaling)
from swells.spectra import (Hm0, bandwidth_nu, cos2s_spreading, jonswap,
                            moment, peak_frequency, pierson_moskowitz,
                            storm_spectrum)

F = np.linspace(0.005, 1.0, 20000)


def test_benchmark_1_generation():
    """U10 = 25 m/s, F = 600 km. Checked by hand in the review."""
    s = sea_state(25.0, 600e3)
    assert s.x_tilde == pytest.approx(9414.384, abs=0.01)
    assert s.f_p == pytest.approx(0.067034, abs=2e-6)
    assert s.T_p == pytest.approx(14.918, abs=0.002)
    assert s.alpha == pytest.approx(0.010153, abs=2e-6)
    assert s.regime == "fetch-limited"


def test_calibrated_spectrum_matches_the_energy_law():
    """m_0 from integrating the spectrum vs g Hs/U^2 = 0.0016 sqrt(x~).

    With calibration on these agree by construction, so this is really a check
    that the internal reference grid used for calibration is wide enough to
    capture the whole spectrum. If it were truncated, the scale factor would
    come out too large and this would fail.
    """
    for U10, fetch in [(15.0, 300e3), (25.0, 600e3), (20.0, 100e3)]:
        S, s = storm_spectrum(F, U10, fetch)
        assert Hm0(F, S) == pytest.approx(
            significant_height_scaling(U10, fetch), rel=0.01)


def test_the_two_jonswap_fits_disagree_and_we_know_by_how_much():
    """Pin down the inconsistency rather than papering over it.

    The alpha and f_p fits imply H_s ~ x~^0.55; the energy fit says x~^0.50.
    So the raw spectrum runs high at large fetch and low at small fetch, and
    the two cross somewhere in between. If a future edit silently "fixes" one
    of the fits, this test notices.
    """
    ratios = {}
    for U10, fetch in [(20.0, 20e3), (25.0, 600e3), (25.0, 5000e3)]:
        S_raw, s = storm_spectrum(F, U10, fetch, calibrate=False)
        ratios[s.x_tilde] = Hm0(F, S_raw) / significant_height_scaling(U10, fetch)

    x_small, x_mid, x_big = sorted(ratios)
    assert ratios[x_small] < ratios[x_mid] < ratios[x_big]
    assert ratios[x_mid] == pytest.approx(1.24, rel=0.02)   # the 24% at x~ ~ 1e4

    # The two laws cross at x~ = (1.60/1.26)^(1/0.05) ~ 1e2, at the very bottom
    # of the JONSWAP validity range. So over the whole range anyone actually
    # uses, the raw spectrum runs high -- never low.
    assert all(r > 1.0 for r in ratios.values())

    # And the ratio must scale as x~^0.05, the difference of the two exponents.
    # This is the sharpest available evidence that the diagnosis is right and
    # not just a coincidence of coefficients.
    predicted = (x_big / x_small) ** 0.05
    assert ratios[x_big] / ratios[x_small] == pytest.approx(predicted, rel=0.005)


def test_benchmark_1_significant_height():
    S, _ = storm_spectrum(F, 25.0, 600e3)
    assert Hm0(F, S) == pytest.approx(9.894, rel=0.02)
    assert moment(F, S, 0) == pytest.approx(6.118, rel=0.04)


def test_gamma_one_recovers_pierson_moskowitz():
    s = sea_state(20.0, 500e3)
    assert np.allclose(jonswap(F, s.alpha, s.f_p, gamma=1.0),
                       pierson_moskowitz(F, f_p=s.f_p, alpha=s.alpha))


def test_peak_of_jonswap_sits_at_fp():
    s = sea_state(25.0, 600e3)
    S = jonswap(F, s.alpha, s.f_p)
    assert peak_frequency(F, S) == pytest.approx(s.f_p, rel=0.02)


def test_growth_saturates_at_the_fully_developed_state():
    """Infinite fetch must not give infinite waves."""
    s = sea_state(20.0, 1e9)
    assert s.regime == "fully-developed"
    assert s.x_tilde == pytest.approx(FULLY_DEVELOPED_FETCH)
    # and the peak frequency should land on the PM value 0.13 g/U10
    assert s.f_p == pytest.approx(0.13 * G / 20.0, rel=0.01)


def test_pierson_moskowitz_fully_developed_height():
    """H_s ~ 0.0246 U10^2 for a fully developed sea."""
    for U10 in [10.0, 20.0, 30.0]:
        S = pierson_moskowitz(F, U10=U10)
        assert Hm0(F, S) == pytest.approx(0.0246 * U10**2, rel=0.08)


def test_duration_can_be_the_binding_limit():
    """A short blow cannot exploit a long fetch."""
    s_long = sea_state(25.0, 2000e3, duration=np.inf)
    s_short = sea_state(25.0, 2000e3, duration=6 * 3600.0)
    assert s_short.regime == "duration-limited"
    assert s_short.f_p > s_long.f_p          # shorter blow, shorter waves
    assert duration_limited_fetch(6 * 3600.0, 25.0) < 2000e3


def test_duration_relation_is_self_consistent():
    """t~ = 68.8 x~^0.67 inverted and re-applied must give back x~."""
    U10 = 20.0
    for x in [1e3, 1e4, 5e4]:
        t = 68.8 * x**0.67 * U10 / G
        assert dimensionless_fetch(duration_limited_fetch(t, U10), U10) == \
            pytest.approx(x, rel=1e-9)


def test_swell_band_is_far_narrower_than_any_wind_sea():
    """The comparison nu is actually good for."""
    narrow = np.exp(-((F - 0.07) ** 2) / (2 * 0.002**2))
    assert bandwidth_nu(F, narrow) == pytest.approx(0.002 / 0.07, rel=1e-3)
    for fetch in [50e3, 500e3, 2000e3]:
        assert bandwidth_nu(F, narrow) < 0.2 * bandwidth_nu(
            F, storm_spectrum(F, 25.0, fetch)[0])


def test_nu_is_sensitive_to_where_you_truncate_the_tail():
    """A defect of nu, documented so nobody trusts it too far.

    An f^-5 spectrum has an m_2 integrand going as f^-3: convergent, but only
    just. Doubling the cutoff still moves nu appreciably. Compare kappa, below,
    which is weighted by S rather than f^2 S and barely moves at all -- which
    is why the grouping calculations use kappa.
    """
    from swells.groups import envelope_correlation

    S, s = storm_spectrum(F, 25.0, 600e3)
    cuts = (3, 6, 12, 25)
    nus = [bandwidth_nu(F, S, f_cut=c * s.f_p) for c in cuts]
    kappas = [envelope_correlation(F[F <= c * s.f_p], S[F <= c * s.f_p])
              for c in cuts]

    # nu drifts by 30% over this range and has still not settled at 25 f_p.
    assert nus[0] < nus[1] < nus[2] < nus[3]
    assert nus[-1] / nus[0] == pytest.approx(1.30, rel=0.05)

    # kappa moves by less than one percent and is converged by 6 f_p.
    assert max(kappas) - min(kappas) < 0.01
    assert kappas[-1] == pytest.approx(kappas[1], rel=0.002)


def test_nu_tracks_fmax_over_fp_not_wave_age():
    """The trap: on a fixed grid an *older* sea reports a larger nu, because
    its peak sits lower so more octaves of tail fall inside the grid. Cut at a
    fixed multiple of f_p instead and the spurious ordering disappears."""
    young, s_y = storm_spectrum(F, 25.0, 50e3)
    old, s_o = storm_spectrum(F, 25.0, 2000e3)
    assert s_o.f_p < s_y.f_p
    assert bandwidth_nu(F, old) > bandwidth_nu(F, young)          # the artefact
    assert bandwidth_nu(F, old, f_cut=4 * s_o.f_p) == pytest.approx(
        bandwidth_nu(F, young, f_cut=4 * s_y.f_p), rel=0.05)      # self-similar


def test_directional_spreading_normalises():
    th = np.linspace(-np.pi, np.pi, 4001)
    for s in [2.0, 10.0, 75.0]:
        D = cos2s_spreading(th, 0.0, s)
        assert trapz(D, th) == pytest.approx(1.0, rel=1e-4)
