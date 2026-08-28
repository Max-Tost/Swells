"""Thin adapter between the browser UI and the swells package.

Runs inside Pyodide. Its only job is to call swells.simulate.run and flatten
the result into something JSON can carry, downsampled to what a canvas can
usefully draw. No physics lives here -- if you find yourself computing
something in this file, it belongs in the package instead.
"""

import json

import numpy as np

from swells import Coast, Storm, propagate, run, surf_report
from swells.util import trapz


def _thin(a, n=600):
    """Downsample an array to at most n points, keeping the endpoints."""
    a = np.asarray(a, dtype=np.float64)
    if a.size <= n:
        return a.tolist()
    idx = np.linspace(0, a.size - 1, n).round().astype(int)
    return a[idx].tolist()


def _energy_band(f, S, lo=0.004, hi=0.996):
    c = np.cumsum(S)
    if c[-1] <= 0:
        return float(f[0]), float(f[-1])
    c = c / c[-1]
    return float(np.interp(lo, c, f)), float(np.interp(hi, c, f))


def simulate(U10, fetch_km, duration_h, distance_km,
             slope, approach_deg, record_hours=1.0, seed=0):
    storm = Storm(U10=U10, fetch_km=fetch_km, duration_h=duration_h,
                  distance_km=distance_km)
    coast = Coast(slope=slope, approach_deg=approach_deg)
    r = run(storm=storm, coast=coast, record_hours=record_hours, seed=seed)

    f = r.f
    f_lo, f_hi = _energy_band(f, r.S_source)
    band = (f >= f_lo * 0.6) & (f <= f_hi * 1.3)
    f_b = f[band]

    # Spectrogram, cropped in both axes and coarsened for the canvas.
    t_lo = 0.6 * propagate.travel_time(f_lo, storm.distance) / 3600.0
    t_hi = 1.25 * propagate.travel_time(f_hi, storm.distance) / 3600.0
    tsel = (r.t_hours >= t_lo) & (r.t_hours <= t_hi)
    if tsel.sum() < 4:
        tsel = np.ones_like(r.t_hours, dtype=bool)
    S_ft = r.S_ft[np.ix_(tsel, band)]
    ti = np.linspace(0, S_ft.shape[0] - 1, min(220, S_ft.shape[0])).round().astype(int)
    fi = np.linspace(0, S_ft.shape[1] - 1, min(140, S_ft.shape[1])).round().astype(int)
    grid = S_ft[np.ix_(ti, fi)]

    # A window of the record short enough to see individual waves in.
    t_rec = r.record["t"]
    win = t_rec <= min(t_rec[-1], 25 * 60.0)

    tab = r.surf["table"]
    ok = np.isfinite(tab["H"])

    b = r.breaking
    out = {
        "report": surf_report(r),
        "storm": {"Tp": r.sea.T_p, "Hs": float(4 * np.sqrt(
            trapz(r.S_source, f))), "regime": r.sea.regime,
            "x_tilde": r.sea.x_tilde},
        "spectrum": {"f": _thin(f_b), "S_source": _thin(r.S_source[band]),
                     "S_peak": _thin(r.S_ft[int(np.argmax(r.Hm0_t))][band])},
        "spectrogram": {"t": r.t_hours[tsel][ti].tolist() if tsel.sum() >= 4
                        else r.t_hours[ti].tolist(),
                        "f": f_b[fi].tolist(),
                        "z": grid.tolist(), "zmax": float(grid.max())},
        "timeline": {"t": _thin(r.t_hours[tsel], 500),
                     "Hm0": _thin(r.Hm0_t[tsel], 500),
                     "Tp": _thin(r.Tp_t[tsel], 500)},
        "record": {"t": _thin(t_rec[win] / 60.0, 2400),
                   "eta": _thin(r.record["eta"][win], 2400),
                   "env": _thin(r.record["envelope"][win], 2400),
                   "Hm0": r.peak["Hm0"]},
        "buoy": {"Hm0": r.peak["Hm0"], "Tp": r.peak["Tp"],
                 "nu": r.peak["nu"], "kappa": r.grouping["kappa"],
                 "waves_per_set": r.grouping["mean_run_length"],
                 "waves_per_set_counted": r.measured.get("mean_run_length", 0.0),
                 "set_interval_min": r.grouping["set_interval"] / 60.0,
                 "peak_day": r.peak["time_h"] / 24.0,
                 "chirp": r.peak["chirp_hz_per_day"],
                 "BFI": r.surf["BFI"],
                 "H_max_hour": r.record["H_hundredth"]},
        "surf": {"x": _thin(tab["x"][ok]), "h": _thin(tab["h"][ok]),
                 "H": _thin(tab["H"][ok]),
                 "break": None if b is None else {
                     "x": b.x, "h": b.h, "H": b.H, "kind": b.kind,
                     "xi0": b.xi0, "gamma_b": b.gamma_b, "peel": b.peel_deg}},
    }
    return json.dumps(out)
