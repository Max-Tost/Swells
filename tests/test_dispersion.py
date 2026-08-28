"""Analytic limits and known values for linear wave theory."""

import numpy as np
import pytest

from swells.constants import G
from swells.dispersion import (deep_group_speed, deep_phase_speed,
                               deep_wavelength, energy_density, group_speed,
                               n_ratio, phase_speed, wavenumber)


def test_dispersion_relation_is_satisfied():
    omega = 2 * np.pi / np.array([4.0, 8.0, 15.0, 20.0])
    for h in [1.0, 5.0, 50.0, 500.0, 4000.0]:
        k = wavenumber(omega, h)
        assert np.allclose(G * k * np.tanh(k * h), omega**2, rtol=1e-10)


def test_deep_water_limit():
    """h -> infinity: k -> omega^2/g, c -> g/omega, c_g -> c/2."""
    omega = 2 * np.pi / 15.0
    h = 10000.0
    assert wavenumber(omega, h) == pytest.approx(omega**2 / G, rel=1e-9)
    assert phase_speed(omega, h) == pytest.approx(G / omega, rel=1e-9)
    assert n_ratio(omega, h) == pytest.approx(0.5, abs=1e-9)
    assert group_speed(omega, h) == pytest.approx(0.5 * G / omega, rel=1e-9)


def test_shallow_water_limit():
    """kh -> 0: c -> sqrt(gh) and c_g -> c, i.e. non-dispersive."""
    T, h = 300.0, 5.0                      # a tsunami-ish period, very shallow
    omega = 2 * np.pi / T
    assert phase_speed(omega, h) == pytest.approx(np.sqrt(G * h), rel=1e-3)
    assert n_ratio(omega, h) == pytest.approx(1.0, abs=1e-3)


def test_deep_water_shortcuts_agree_with_general_formulas():
    f = np.array([0.05, 0.08, 0.12])
    omega = 2 * np.pi * f
    h = 6000.0
    assert np.allclose(deep_phase_speed(f), phase_speed(omega, h), rtol=1e-8)
    assert np.allclose(deep_group_speed(f), group_speed(omega, h), rtol=1e-8)
    assert np.allclose(deep_wavelength(f), 2 * np.pi / wavenumber(omega, h), rtol=1e-8)


def test_hand_checkable_anchors():
    """A 20 s swell runs at 56 km/h, so 10,000 km takes 7.4 days."""
    cg = deep_group_speed(1 / 20.0)
    assert cg == pytest.approx(15.61, abs=0.01)
    assert cg * 3.6 == pytest.approx(56.2, abs=0.1)
    days = 1.0e7 / cg / 86400.0
    assert days == pytest.approx(7.41, abs=0.02)

    # A 15 s wave is 351 m long in deep water.
    assert deep_wavelength(1 / 15.0) == pytest.approx(351.2, abs=0.2)


def test_energy_density():
    """E = rho g H^2 / 8. A 2 m wave carries about 5 kJ/m^2."""
    assert energy_density(2.0) == pytest.approx(1025 * G * 4.0 / 8.0, rel=1e-12)
    assert energy_density(2.0) / 1e3 == pytest.approx(5.03, abs=0.02)


def test_benchmark_3_deep_water_baseline():
    """From the research draft, T = 15 s."""
    T = 15.0
    omega = 2 * np.pi / T
    assert omega == pytest.approx(0.418879, abs=1e-6)
    assert omega**2 / G == pytest.approx(0.0178919, abs=1e-7)
    assert deep_wavelength(1 / T) == pytest.approx(351.17, abs=0.05)
    assert deep_phase_speed(1 / T) == pytest.approx(23.411, abs=0.01)
    assert deep_group_speed(1 / T) == pytest.approx(11.706, abs=0.01)


def test_benchmark_3_intermediate_depth_row():
    """The h = 10 m row of the transformation table."""
    omega = 2 * np.pi / 15.0
    h = 10.0
    k = wavenumber(omega, h)
    assert k == pytest.approx(0.04360, abs=1e-5)
    assert 2 * np.pi / k == pytest.approx(144.10, abs=0.05)
    assert phase_speed(omega, h) == pytest.approx(9.607, abs=0.005)
    assert group_speed(omega, h) == pytest.approx(9.049, abs=0.005)
