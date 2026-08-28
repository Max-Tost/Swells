"""Generate every figure used in the lessons.

Run from the repository root:

    .venv/bin/python figures/make_figures.py

Nothing in lessons/ is drawn by hand. If a plot disagrees with the text, one of
them is wrong and this script is the arbiter.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from swells import Coast, Storm, groups, nearshore, propagate, run, spectra, synth
from swells.constants import G
from swells.dispersion import (deep_group_speed, deep_wavelength, group_speed,
                               n_ratio, phase_speed, wavenumber)
from swells.util import trapz

OUT = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
})
INK = "#1b3a5c"
ACCENT = "#c0492b"
MUTED = "#7f8c9a"


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print("  ", name)


# --- 01/02: dispersion and group velocity ---------------------------------

def fig_dispersion():
    kh = np.logspace(-1.5, 1.2, 500)
    h = 100.0
    k = kh / h
    omega = np.sqrt(G * k * np.tanh(kh))
    c = omega / k
    n = 0.5 * (1 + 2 * kh / np.sinh(2 * kh))

    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.2))
    ax[0].semilogx(kh, c / np.sqrt(G * h), color=INK, label="$c/\\sqrt{gh}$")
    ax[0].semilogx(kh, n * c / np.sqrt(G * h), color=ACCENT, label="$c_g/\\sqrt{gh}$")
    ax[0].axvline(np.pi / 10, color=MUTED, ls=":", lw=1)
    ax[0].axvline(np.pi, color=MUTED, ls=":", lw=1)
    ax[0].text(np.pi / 30, 0.9, "shallow", color=MUTED, fontsize=8)
    ax[0].text(6, 0.9, "deep", color=MUTED, fontsize=8)
    ax[0].set_xlabel("$kh$"); ax[0].set_ylabel("speed / $\\sqrt{gh}$")
    ax[0].legend(frameon=False)

    ax[1].semilogx(kh, n, color=INK)
    ax[1].axhline(0.5, color=MUTED, ls="--", lw=1)
    ax[1].axhline(1.0, color=MUTED, ls="--", lw=1)
    ax[1].set_xlabel("$kh$"); ax[1].set_ylabel("$n = c_g/c$")
    ax[1].set_title("energy travels at half the crest speed in deep water",
                    fontsize=9, color=MUTED)
    save(fig, "fig01_dispersion.png")


def fig_group_conveyor():
    """Crests overtaking their own group: c_g = c/2 made visible."""
    f0, df = 0.08, 0.004
    x = np.linspace(0, 2500, 3000)
    fig, ax = plt.subplots(figsize=(8.4, 3.6))
    for i, t in enumerate([0, 40, 80]):
        k1 = (2 * np.pi * (f0 + df)) ** 2 / G
        k2 = (2 * np.pi * (f0 - df)) ** 2 / G
        w1, w2 = 2 * np.pi * (f0 + df), 2 * np.pi * (f0 - df)
        eta = np.cos(k1 * x - w1 * t) + np.cos(k2 * x - w2 * t)
        env = 2 * np.abs(np.cos(0.5 * ((k1 - k2) * x - (w1 - w2) * t)))
        off = -i * 5
        ax.plot(x, eta + off, color=INK, lw=0.7)
        ax.plot(x, env + off, color=ACCENT, lw=1.3)
        ax.plot(x, -env + off, color=ACCENT, lw=1.3)
        ax.text(20, off + 2.6, f"$t = {t}$ s", fontsize=8, color=MUTED)

    cg = deep_group_speed(f0)
    c = 2 * cg
    for t, m in [(0, "o"), (40, "o"), (80, "o")]:
        ax.plot(600 + cg * t, -[0, 40, 80].index(t) * 5, m, color=ACCENT, ms=5)
        ax.plot(600 + c * t, -[0, 40, 80].index(t) * 5, "x", color=INK, ms=6)
    ax.set_xlabel("x (m)"); ax.set_yticks([])
    ax.set_title("red = envelope (moves at $c_g$), black x = one crest "
                 "(moves at $2c_g$)", fontsize=9, color=MUTED)
    save(fig, "fig02_group_conveyor.png")


# --- 04: spectra -----------------------------------------------------------

def fig_spectra():
    f = np.linspace(0.02, 0.5, 2000)
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.2))

    for fetch, col in [(50e3, "#9ec5e8"), (200e3, "#5b93c7"), (600e3, INK),
                       (5000e3, ACCENT)]:
        S, s = spectra.storm_spectrum(f, 25.0, fetch)
        ax[0].plot(f, S, color=col,
                   label=f"{fetch/1e3:.0f} km  ($T_p$={s.T_p:.0f} s)")
    ax[0].set_xlabel("f (Hz)"); ax[0].set_ylabel("S(f)  (m$^2$/Hz)")
    ax[0].set_title("growth with fetch, $U_{10}=25$ m/s", fontsize=9, color=MUTED)
    ax[0].legend(frameon=False, fontsize=7.5)

    s = spectra.sea_state(25.0, 600e3)
    for gamma, col in [(1.0, MUTED), (3.3, INK), (7.0, ACCENT)]:
        ax[1].plot(f, spectra.jonswap(f, s.alpha, s.f_p, gamma=gamma),
                   color=col, label=f"$\\gamma$ = {gamma}")
    ax[1].set_xlim(0.03, 0.2)
    ax[1].set_xlabel("f (Hz)")
    ax[1].set_title("peak enhancement ($\\gamma$=1 is Pierson-Moskowitz)",
                    fontsize=9, color=MUTED)
    ax[1].legend(frameon=False)
    save(fig, "fig04_spectra.png")


def fig_jonswap_inconsistency():
    U10 = 25.0
    x = np.logspace(2, 4.4, 60)
    fetch = x * U10**2 / G
    f = np.linspace(0.01, 1.2, 12000)
    raw, law = [], []
    for F_ in fetch:
        S, s = spectra.storm_spectrum(f, U10, F_, calibrate=False)
        raw.append(spectra.Hm0(f, S) * G / U10**2)
        law.append(0.0016 * np.sqrt(s.x_tilde))

    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ax.loglog(x, raw, color=ACCENT, label="from the $\\alpha$, $f_p$ fits  ($\\propto \\tilde x^{0.55}$)")
    ax.loglog(x, law, color=INK, label="from the energy fit  ($\\propto \\tilde x^{0.50}$)")
    ax.set_xlabel("dimensionless fetch  $\\tilde x$")
    ax.set_ylabel("$g H_s / U_{10}^2$")
    ax.set_title("two JONSWAP results that disagree", fontsize=9, color=MUTED)
    ax.legend(frameon=False, fontsize=8)
    save(fig, "fig04_inconsistency.png")


# --- 05/06: the crossing ---------------------------------------------------

def energy_band(f, S, lo=0.005, hi=0.995):
    """The frequency range holding the bulk of the energy.

    Needed because a JONSWAP tail runs to 1 Hz with nothing in it, and 1 Hz
    takes four months to cross an ocean. Plotting the full grid gives a mostly
    empty picture.
    """
    c = np.cumsum(S) / np.sum(S)
    return float(np.interp(lo, c, f)), float(np.interp(hi, c, f))


def fig_spectrogram():
    D = 8.0e6
    f = np.linspace(0.02, 0.25, 700)
    S0, sea = spectra.storm_spectrum(f, 24.0, 900e3, 40 * 3600.0)
    f_lo, f_hi = energy_band(f, S0)
    t = np.linspace(0.6 * propagate.travel_time(f_lo, D),
                    1.25 * propagate.travel_time(f_hi, D), 900)
    S_ft = propagate.spectrogram(f, S0, D, 40 * 3600.0, 900e3, t)
    Hm0_t, _ = propagate.sea_state_at(f, S_ft)
    ridge = propagate.spectrogram_ridge(f, S_ft)
    df_dt, t0 = propagate.fit_chirp(t, ridge, weights=Hm0_t**2)
    D_fit = propagate.distance_from_chirp(df_dt)

    fig, ax = plt.subplots(1, 2, figsize=(8.8, 3.4),
                           gridspec_kw={"width_ratios": [2, 1]})
    ax[0].pcolormesh(t / 86400, f, S_ft.T, shading="auto", cmap="magma")
    ax[0].plot(t / 86400, df_dt * t + (-df_dt * t0), color="white", lw=1.2, ls="--")
    ax[0].set_ylim(f_lo * 0.7, f_hi * 1.15)
    ax[0].set_xlabel("days since the storm"); ax[0].set_ylabel("f (Hz)")
    ax[0].set_title(f"true D = {D/1e3:.0f} km,  fitted D = {D_fit/1e3:.0f} km",
                    fontsize=9, color=MUTED)

    ax[1].plot(t / 86400, Hm0_t, color=INK)
    ax[1].set_xlabel("days"); ax[1].set_ylabel("$H_{m0}$ (m)")
    ax[1].set_title("long period first", fontsize=9, color=MUTED)
    save(fig, "fig05_spectrogram.png")


def fig_attenuation():
    f = np.linspace(0.02, 0.25, 600)
    S0, _ = spectra.storm_spectrum(f, 25.0, 600e3, 36 * 3600.0)
    D = np.logspace(5.7, 7.2, 40)
    H, band = [], []
    for d in D:
        t = np.linspace(0.4 * propagate.travel_time(f.min(), d),
                        1.2 * propagate.travel_time(f.max(), d), 500)
        S_ft = propagate.spectrogram(f, S0, d, 36 * 3600.0, 600e3, t)
        H.append(np.max(propagate.sea_state_at(f, S_ft)[0]))
        band.append(propagate.arrival_bandwidth(0.07, d, 36 * 3600.0, 600e3))

    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.2))
    ax[0].loglog(D / 1e3, H, color=INK)
    ax[0].set_xlabel("distance (km)"); ax[0].set_ylabel("peak $H_{m0}$ (m)")
    ax[0].set_title("smaller with distance", fontsize=9, color=MUTED)
    ax[1].loglog(D / 1e3, band, color=ACCENT)
    ax[1].set_xlabel("distance (km)"); ax[1].set_ylabel("$\\Delta f$ (Hz)")
    ax[1].set_title("and cleaner", fontsize=9, color=MUTED)
    save(fig, "fig06_attenuation.png")


# --- 07: sets --------------------------------------------------------------

def fig_sets():
    r = run(record_hours=1.2, dt=0.4, seed=3)
    t = r.record["t"] / 60.0
    eta, env = r.record["eta"], r.record["envelope"]
    Hm0 = r.peak["Hm0"]

    fig, ax = plt.subplots(2, 1, figsize=(8.8, 5.0))
    m = t < 25
    ax[0].plot(t[m], eta[m], color=INK, lw=0.6)
    ax[0].plot(t[m], env[m], color=ACCENT, lw=1.2)
    ax[0].plot(t[m], -env[m], color=ACCENT, lw=1.2)
    ax[0].axhline(Hm0 / 2, color=MUTED, ls=":", lw=1)
    ax[0].axhline(-Hm0 / 2, color=MUTED, ls=":", lw=1)
    ax[0].set_xlabel("minutes"); ax[0].set_ylabel("$\\eta$ (m)")
    ax[0].set_title(
        f"$T_p$ = {r.peak['Tp']:.0f} s,  predicted {r.grouping['mean_run_length']:.1f} "
        f"waves per set,  counted {r.measured['mean_run_length']:.1f}",
        fontsize=9, color=MUTED)

    # The upper panel is 25 minutes because that is what you can see. The
    # histogram needs far more waves than that, so count over 60 hours.
    t_long, eta_long = synth.surface_elevation(
        r.f, r.S_ft[np.argmax(r.Hm0_t)], 60 * 3600.0, 0.5, seed=11)
    long_stats = groups.measured_group_statistics(t_long, eta_long)

    j, P = groups.run_length_distribution(r.grouping["p22"], 10)
    runs = long_stats["runs"]
    counted = np.array([(runs == jj).sum() for jj in j], float)
    counted /= max(counted.sum(), 1)
    w = 0.38
    ax[1].bar(j - w / 2, P, w, color=INK, label="predicted (bivariate Rayleigh)")
    ax[1].bar(j + w / 2, counted, w, color=ACCENT, label="counted in the record")
    ax[1].set_xlabel("waves in a run"); ax[1].set_ylabel("probability")
    ax[1].legend(frameon=False)
    ax[1].set_title(f"{long_stats['n_groups']} runs counted over 60 hours",
                    fontsize=9, color=MUTED)
    save(fig, "fig07_sets.png")


def fig_grouping_vs_bandwidth():
    F = np.linspace(0.02, 0.4, 1500)
    widths = np.linspace(0.002, 0.03, 22)
    kappa, jbar = [], []
    for w in widths:
        S = np.exp(-((F - 0.07) ** 2) / (2 * w**2))
        S = S / trapz(S, F)
        g = groups.group_statistics(F, S)
        kappa.append(g["kappa"]); jbar.append(g["mean_run_length"])

    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ax.plot(widths / 0.07, jbar, color=INK, marker="o", ms=3)
    ax.axhline(1.156, color=MUTED, ls="--", lw=1)
    ax.text(0.28, 1.20, "uncorrelated waves: 1.16", color=MUTED, fontsize=8)
    ax.set_xlabel("relative bandwidth  $\\sigma_f / f_p$")
    ax.set_ylabel("mean waves per set")
    ax.set_title("narrow spectra are groupy", fontsize=9, color=MUTED)
    save(fig, "fig07_grouping.png")


# --- 09/10: nearshore ------------------------------------------------------

def fig_shoaling():
    T = 15.0
    omega = 2 * np.pi / T
    h = np.linspace(200, 1.0, 800)
    Ks = nearshore.shoaling_coefficient(omega, h)
    th0 = np.radians(35.0)
    th = nearshore.snell_angle(omega, h, th0)
    Kr = nearshore.refraction_coefficient(th0, th)

    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.2))
    ax[0].plot(h, Ks, color=INK, label="$K_s$")
    ax[0].plot(h, Kr, color=ACCENT, label="$K_r$ (35$^\\circ$ approach)")
    ax[0].plot(h, Ks * Kr, color="black", lw=1.4, ls="--", label="$K_sK_r$")
    ax[0].axhline(1, color=MUTED, lw=0.8)
    ax[0].set_xscale("log"); ax[0].invert_xaxis()
    ax[0].set_xlabel("depth (m)"); ax[0].set_ylabel("coefficient")
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].set_title("$K_s$ dips to 0.913 before Green's law", fontsize=9, color=MUTED)

    ax[1].plot(h, np.degrees(th), color=INK)
    ax[1].set_xscale("log"); ax[1].invert_xaxis()
    ax[1].set_xlabel("depth (m)"); ax[1].set_ylabel("angle to the contours (deg)")
    ax[1].set_title("everything turns to face the beach", fontsize=9, color=MUTED)
    save(fig, "fig09_shoaling.png")


def fig_breaking():
    fig, ax = plt.subplots(1, 2, figsize=(8.8, 3.4))

    for slope, col, lab in [(0.01, "#9ec5e8", "1:100"), (0.02, "#5b93c7", "1:50"),
                            (0.05, INK, "1:20"), (0.10, ACCENT, "1:10")]:
        x, h = nearshore.plane_beach(slope, 3000.0, 900)
        tab, ev = nearshore.transform_transect(15.0, 2.0, 0.0, x, h, tan_beta=slope)
        ok = np.isfinite(tab["H"])
        ax[0].plot(tab["h"][ok], tab["H"][ok], color=col,
                   label=f"{lab}  $\\xi_0$={ev.xi0:.2f}  {ev.kind}")
        ax[0].plot(ev.h, ev.H, "o", color=col, ms=5)
    ax[0].set_xscale("log"); ax[0].invert_xaxis()
    ax[0].set_xlim(60, 0.5)
    ax[0].set_xlabel("depth (m)"); ax[0].set_ylabel("H (m)")
    ax[0].legend(frameon=False, fontsize=7.5)
    ax[0].set_title("2 m of 15 s swell on four beaches", fontsize=9, color=MUTED)

    xi = np.logspace(-1, 1, 400)
    for lo, hi, name, col in [(0, 0.5, "spilling", "#9ec5e8"),
                              (0.5, 3.3, "plunging", "#5b93c7"),
                              (3.3, 3.8, "collapsing", INK),
                              (3.8, 10, "surging", ACCENT)]:
        ax[1].axvspan(max(lo, 0.1), hi, color=col, alpha=0.35)
        ax[1].text(np.sqrt(max(lo, 0.1) * hi), 0.5, name, rotation=90,
                   ha="center", va="center", fontsize=8)
    ax[1].set_xscale("log"); ax[1].set_xlim(0.1, 10); ax[1].set_yticks([])
    ax[1].set_xlabel("Iribarren number  $\\xi_0$")
    ax[1].set_title("one number decides how it falls over", fontsize=9, color=MUTED)
    save(fig, "fig10_breaking.png")


# --- 11: the chain ---------------------------------------------------------

def fig_chain():
    r = run(record_hours=0.8, seed=5)
    fig, ax = plt.subplots(2, 2, figsize=(9.0, 5.6))

    ax[0, 0].plot(r.f, r.S_source, color=INK)
    ax[0, 0].set_xlim(0, 0.25)
    ax[0, 0].set_xlabel("f (Hz)"); ax[0, 0].set_ylabel("S (m$^2$/Hz)")
    ax[0, 0].set_title("1. in the storm", fontsize=9, color=MUTED)

    f_lo, f_hi = energy_band(r.f, r.S_source)
    ax[0, 1].pcolormesh(r.t_hours / 24, r.f, r.S_ft.T, shading="auto", cmap="magma")
    ax[0, 1].set_ylim(f_lo * 0.7, f_hi * 1.15)
    ax[0, 1].set_xlim(0.6 * propagate.travel_time(f_lo, r.storm.distance) / 86400,
                      1.25 * propagate.travel_time(f_hi, r.storm.distance) / 86400)
    ax[0, 1].set_xlabel("days"); ax[0, 1].set_ylabel("f (Hz)")
    ax[0, 1].set_title("2. the crossing", fontsize=9, color=MUTED)

    t = r.record["t"] / 60
    m = t < 15
    ax[1, 0].plot(t[m], r.record["eta"][m], color=INK, lw=0.6)
    ax[1, 0].plot(t[m], r.record["envelope"][m], color=ACCENT, lw=1.1)
    ax[1, 0].set_xlabel("minutes"); ax[1, 0].set_ylabel("$\\eta$ (m)")
    ax[1, 0].set_title("3. at the buoy: sets", fontsize=9, color=MUTED)

    tab = r.surf["table"]
    ok = np.isfinite(tab["H"])
    ax[1, 1].plot(tab["x"][ok], tab["H"][ok], color=INK)
    ax[1, 1].plot(tab["x"][ok], -tab["h"][ok] / 3, color=MUTED, lw=1)
    if r.breaking:
        ax[1, 1].axvline(r.breaking.x, color=ACCENT, ls="--", lw=1)
        ax[1, 1].text(r.breaking.x, 0.2,
                      f" {r.breaking.kind}\n {r.breaking.H:.1f} m", fontsize=8,
                      color=ACCENT)
    ax[1, 1].invert_xaxis()
    ax[1, 1].set_xlabel("distance offshore (m)"); ax[1, 1].set_ylabel("H (m) / depth/3")
    ax[1, 1].set_title("4. the surf zone", fontsize=9, color=MUTED)
    save(fig, "fig11_chain.png")


if __name__ == "__main__":
    print("writing figures to", OUT)
    for fn in [fig_dispersion, fig_group_conveyor, fig_spectra,
               fig_jonswap_inconsistency, fig_spectrogram, fig_attenuation,
               fig_sets, fig_grouping_vs_bandwidth, fig_shoaling, fig_breaking,
               fig_chain]:
        fn()
    print("done")
