"""Shoaling, refraction, breaking -- and conservation of energy flux."""

import numpy as np
import pytest

from swells.constants import G
from swells.dispersion import deep_wavelength, energy_flux, group_speed
from swells.nearshore import (breaker_type, goda_gamma, iribarren, miche_limit,
                              plane_beach, refraction_coefficient, set_down,
                              set_up_slope, shoaling_coefficient, snell_angle,
                              transform_transect)


def test_shoaling_conserves_energy_flux():
    """The defining property: E c_g is the same at every depth."""
    T = 15.0
    omega = 2 * np.pi / T
    H0 = 2.0
    F0 = energy_flux(H0, omega, 5000.0)
    for h in [200.0, 80.0, 30.0, 12.0, 6.0, 3.0]:
        H = H0 * shoaling_coefficient(omega, h)
        assert energy_flux(H, omega, h) == pytest.approx(F0, rel=1e-6)


def test_shoaling_coefficient_dips_below_one():
    """K_s has a minimum of about 0.913 before Green's law takes over."""
    omega = 2 * np.pi / 15.0
    h = np.linspace(3.0, 300.0, 2000)
    Ks = shoaling_coefficient(omega, h)
    assert Ks.min() == pytest.approx(0.913, abs=0.01)


def test_greens_law_in_shallow_water():
    """K_s -> h^(-1/4)."""
    omega = 2 * np.pi / 60.0            # long period, so 2 m really is shallow
    r = shoaling_coefficient(omega, 1.0) / shoaling_coefficient(omega, 16.0)
    assert r == pytest.approx(16.0**0.25, rel=0.02)


def test_snell_turns_waves_toward_the_beach():
    omega = 2 * np.pi / 15.0
    th0 = np.radians(45.0)
    angles = [np.degrees(snell_angle(omega, h, th0)) for h in [200, 50, 20, 5, 1]]
    assert all(b < a for a, b in zip(angles, angles[1:]))
    assert angles[-1] < 10.0


def test_refraction_coefficient_reduces_height_for_oblique_waves():
    omega = 2 * np.pi / 15.0
    th0 = np.radians(40.0)
    th = snell_angle(omega, 10.0, th0)
    assert refraction_coefficient(th0, th) < 1.0


def test_benchmark_3_table_row_h10():
    """Reproduce the h = 10 m row from the research draft."""
    T, H0, th0_deg = 15.0, 3.0, 30.0
    omega = 2 * np.pi / T
    h = 10.0
    th0 = np.radians(th0_deg)
    th = snell_angle(omega, h, th0)
    Ks = shoaling_coefficient(omega, h)
    Kr = refraction_coefficient(th0, th)
    assert np.degrees(th) == pytest.approx(11.83, abs=0.02)
    assert Ks == pytest.approx(1.137, abs=0.003)
    assert Kr == pytest.approx(0.941, abs=0.003)
    assert H0 * Ks * Kr == pytest.approx(3.21, abs=0.01)
    assert miche_limit(omega, h) == pytest.approx(8.40, abs=0.03)


def test_miche_shallow_limit():
    """H_max -> 0.142 * 2 pi * h = 0.892 h."""
    omega = 2 * np.pi / 120.0
    h = 2.0
    assert miche_limit(omega, h) / h == pytest.approx(0.892, rel=0.02)


def test_goda_reduces_to_mccowan_on_a_flat_bed():
    """Flat beach, shallow water: gamma_b should land near 0.78-0.8."""
    L0 = deep_wavelength(1 / 12.0)
    g = goda_gamma(1.5, L0, 0.0)
    assert 0.7 < g < 0.9


def test_goda_gamma_rises_with_slope():
    L0 = deep_wavelength(1 / 12.0)
    assert goda_gamma(1.5, L0, 0.10) > goda_gamma(1.5, L0, 0.01)


def test_iribarren_classification():
    L0 = deep_wavelength(1 / 14.0)
    assert breaker_type(iribarren(2.0, L0, 0.01)) == "spilling"
    assert breaker_type(iribarren(2.0, L0, 0.05)) == "plunging"
    assert breaker_type(iribarren(2.0, L0, 0.30)) == "collapsing"
    assert breaker_type(iribarren(2.0, L0, 0.45)) == "surging"


def test_set_down_is_negative_and_small():
    omega = 2 * np.pi / 15.0
    sd = set_down(2.0, omega, 8.0)
    assert sd < 0
    assert abs(sd) < 0.15


def test_set_up_slope_is_about_a_fifth_of_the_beach_slope():
    assert set_up_slope(0.02) / 0.02 == pytest.approx(0.186, abs=0.005)


def test_transect_breaks_and_is_monotone_inside_the_surf_zone():
    x, h = plane_beach(0.02)
    table, event = transform_transect(15.0, 3.0, 30.0, x, h, tan_beta=0.02)
    assert event is not None
    assert 3.0 < event.h < 6.0
    assert event.kind == "spilling"
    inside = table["h"] < event.h
    Hs = table["H"][inside & np.isfinite(table["H"])]
    assert np.all(np.diff(Hs) <= 1e-9)          # height only decreases once broken


def test_steeper_beach_gives_a_more_violent_breaker():
    x1, h1 = plane_beach(0.01)
    x2, h2 = plane_beach(0.08)
    _, e1 = transform_transect(14.0, 2.0, 0.0, x1, h1, tan_beta=0.01)
    _, e2 = transform_transect(14.0, 2.0, 0.0, x2, h2, tan_beta=0.08)
    assert e2.xi0 > e1.xi0
    assert e1.kind == "spilling" and e2.kind == "plunging"


def test_head_on_swell_closes_out():
    x, h = plane_beach(0.02)
    _, e = transform_transect(15.0, 2.0, 0.0, x, h, tan_beta=0.02)
    assert e.peel_deg == pytest.approx(0.0, abs=1e-9)
