"""Shoaling, refraction, set-up, and breaking.

Derived in lessons/09-feeling-the-bottom.md and lessons/10-breaking.md.

Everything here assumes a mild slope, so that the wave has time to adjust to
the local depth and we can treat it as locally periodic. On a steep reef that
assumption fails and you need a phase-resolving model.
"""

from dataclasses import dataclass

import numpy as np

from .constants import G
from .dispersion import (deep_group_speed, deep_phase_speed, deep_wavelength,
                         group_speed, n_ratio, phase_speed, wavenumber)


# ---------------------------------------------------------------------------
# Non-breaking transformation
# ---------------------------------------------------------------------------

def shoaling_coefficient(omega, h):
    """K_s = sqrt(c_g0 / c_g), from conservation of energy flux E c_g.

    Dips to a minimum of about 0.913 in intermediate depth -- the wave briefly
    gets *smaller* -- then grows like h^(-1/4) (Green's law) in shallow water.
    """
    cg0 = 0.5 * (G / omega)          # deep-water group speed for this omega
    return np.sqrt(cg0 / group_speed(omega, h))


def snell_angle(omega, h, theta0, c0=None):
    """sin(theta)/c = const. Crests swing round to face the beach."""
    if c0 is None:
        c0 = G / omega
    s = np.clip(phase_speed(omega, h) / c0 * np.sin(theta0), -1.0, 1.0)
    return np.arcsin(s)


def refraction_coefficient(theta0, theta):
    """K_r = sqrt(cos theta0 / cos theta), from the widening of the ray tube.

    Blows up at a caustic, where neighbouring rays cross. A 1-D model cannot
    see that coming; see the caveats in lesson 09.
    """
    return np.sqrt(np.cos(theta0) / np.cos(theta))


def set_down(H, omega, h):
    """Mean surface depression under a shoaling wave train.

        eta_bar = -(1/8) H^2 k / sinh(2 k h)

    Radiation stress rises as the wave shoals, and the mean surface tilts down
    to supply the pressure gradient that balances it (Longuet-Higgins &
    Stewart 1964). Small -- centimetres -- but it is the same mechanism that
    produces the much larger set-up inside the surf zone.
    """
    k = wavenumber(omega, h)
    return -(H**2) * k / (8.0 * np.sinh(np.clip(2.0 * k * h, 1e-9, 350.0)))


def set_up_slope(tan_beta, gamma_b=0.78):
    """d(eta_bar)/dx inside the surf zone.

        d eta_bar / dx = tan(beta) / (1 + 8/(3 gamma_b^2))

    For gamma_b = 0.78 the factor is 0.186, so the mean water level climbs at
    about a fifth of the beach slope. On a 1:50 beach with a 2 m breaker that
    is roughly 0.2 m of extra water at the shoreline.
    """
    return tan_beta / (1.0 + 8.0 / (3.0 * gamma_b**2))


def radiation_stress_xx(E, omega, h, theta=0.0):
    """S_xx = E [ n (1 + cos^2 theta) - 1/2 ]."""
    n = n_ratio(omega, h)
    return E * (n * (1.0 + np.cos(theta) ** 2) - 0.5)


# ---------------------------------------------------------------------------
# Breaking
# ---------------------------------------------------------------------------

STOKES_STEEPNESS = 0.142     # (H/L)_max in deep water; crest angle 120 degrees


def miche_limit(omega, h):
    """H_max = 0.142 L tanh(k h) (Miche 1944).

    Interpolates between the deep-water steepness limit and a depth limit. In
    shallow water tanh(kh) -> kh = 2 pi h/L and it collapses to
    H_max = 0.142 * 2 pi * h = 0.892 h.
    """
    k = wavenumber(omega, h)
    L = 2.0 * np.pi / k
    return STOKES_STEEPNESS * L * np.tanh(k * h)


def goda_gamma(h, L0, tan_beta):
    """Goda's (2010) breaker index, which knows about beach slope.

        gamma_b = 0.17 (L0/h) {1 - exp[-1.5 pi (h/L0)(1 + 15 tan(beta)^{4/3})]}

    Steeper beaches let the wave carry further before it collapses, so
    gamma_b rises above the flat-bed value. McCowan's classic 0.78 is the
    tan(beta) -> 0 limit; on a steep beach the observed ratio reaches 1.3.
    """
    h = np.asarray(h, dtype=np.float64)
    arg = -1.5 * np.pi * (h / L0) * (1.0 + 15.0 * np.abs(tan_beta) ** (4.0 / 3.0))
    return 0.17 * (L0 / h) * (1.0 - np.exp(arg))


def iribarren(H0, L0, tan_beta):
    """The surf similarity parameter xi_0 = tan(beta) / sqrt(H0/L0).

    A ratio of two slopes: the beach's, and the wave's own steepness. That one
    number decides how the wave breaks.
    """
    return tan_beta / np.sqrt(H0 / L0)


def breaker_type(xi0):
    """Spilling / plunging / collapsing / surging, by Iribarren number.

    The boundaries are conventional and fuzzy; different authors quote 0.4-0.5
    and 3.0-3.3. Nothing changes discontinuously in the real ocean.
    """
    if xi0 < 0.5:
        return "spilling"
    if xi0 < 3.3:
        return "plunging"
    if xi0 < 3.8:
        return "collapsing"
    return "surging"


def peel_angle(theta_b):
    """Angle between the breaking crest and the shoreline, in degrees.

    The breaker peels along the beach at c_b / sin(peel angle). A small peel
    angle gives a fast, makeable wave; a peel angle of zero is a closeout,
    where the whole crest breaks at once and there is nowhere to go.
    """
    return np.degrees(np.abs(theta_b))


# ---------------------------------------------------------------------------
# The transect
# ---------------------------------------------------------------------------

@dataclass
class Breaking:
    x: float
    h: float
    H: float
    L: float
    theta_deg: float
    gamma_b: float
    xi0: float
    kind: str
    peel_deg: float


def transform_transect(T, H0, theta0_deg, x, h, tan_beta=None):
    """March a wave in over a depth profile h(x), and find where it breaks.

    x should run from offshore to onshore. Returns (table, Breaking or None),
    where table is a dict of arrays.

    Outside the surf zone H = H0 K_s K_r. Once the limit is reached the wave is
    taken to be depth-limited from there in, H = min(Miche, Goda) -- a crude
    but standard stand-in for a real dissipation model such as Battjes &
    Janssen (1978).
    """
    omega = 2.0 * np.pi / T
    theta0 = np.radians(theta0_deg)
    c0 = deep_phase_speed(1.0 / T)
    cg0 = deep_group_speed(1.0 / T)
    L0 = deep_wavelength(1.0 / T)

    x = np.asarray(x, dtype=np.float64)
    h = np.asarray(h, dtype=np.float64)
    if tan_beta is None:
        tan_beta = np.abs(np.gradient(h, x))
    else:
        tan_beta = np.full_like(h, float(tan_beta))

    ok = h > 0.05
    k = np.full_like(h, np.nan)
    k[ok] = wavenumber(omega, h[ok])
    L = 2.0 * np.pi / k
    c = omega / k
    cg = np.full_like(h, np.nan)
    cg[ok] = group_speed(omega, h[ok])

    theta = np.full_like(h, np.nan)
    theta[ok] = snell_angle(omega, h[ok], theta0, c0)
    Ks = np.sqrt(cg0 / cg)
    Kr = np.sqrt(np.cos(theta0) / np.cos(theta))
    H_linear = H0 * Ks * Kr

    H_miche = np.full_like(h, np.nan)
    H_miche[ok] = miche_limit(omega, h[ok])
    H_goda = np.full_like(h, np.nan)
    H_goda[ok] = goda_gamma(h[ok], L0, np.mean(tan_beta[ok])) * h[ok]
    H_limit = np.fmin(H_miche, H_goda)

    exceeds = ok & (H_linear >= H_limit)
    idx = int(np.argmax(exceeds)) if exceeds.any() else None

    H = H_linear.copy()
    event = None
    if idx is not None:
        H[idx:] = np.fmin(H_linear[idx:], H_limit[idx:])
        xi0 = iribarren(H0, L0, float(tan_beta[idx]))
        event = Breaking(
            x=float(x[idx]), h=float(h[idx]), H=float(H_limit[idx]),
            L=float(L[idx]), theta_deg=float(np.degrees(theta[idx])),
            gamma_b=float(H_limit[idx] / h[idx]), xi0=float(xi0),
            kind=breaker_type(xi0), peel_deg=peel_angle(theta[idx]),
        )

    eta_bar = np.full_like(h, np.nan)
    eta_bar[ok] = set_down(H[ok], omega, h[ok])

    table = {"x": x, "h": h, "k": k, "L": L, "c": c, "cg": cg,
             "theta_deg": np.degrees(theta), "Ks": Ks, "Kr": Kr,
             "H_linear": H_linear, "H": H, "H_miche": H_miche,
             "H_goda": H_goda, "set_down": eta_bar, "tan_beta": tan_beta}
    return table, event


def plane_beach(slope, x_offshore=2000.0, n=400):
    """A simple plane beach: h = slope * x, x running from offshore to zero."""
    x = np.linspace(x_offshore, 0.0, n)
    return x, slope * x
