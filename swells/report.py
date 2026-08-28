"""The surf report: the whole simulation in words a person would say out loud."""

import numpy as np


def surf_report(result):
    r = result
    s, c, p = r.storm, r.coast, r.peak
    g, m, b = r.grouping, r.measured, r.breaking

    L = []
    A = L.append

    A("GENERATION")
    A(f"  {s.U10:.0f} m/s wind, {s.fetch_km:.0f} km fetch, {s.duration_h:.0f} h")
    A(f"  regime            {r.sea.regime}")
    A(f"  peak period       {r.sea.T_p:.1f} s  (f_p = {r.sea.f_p:.4f} Hz)")
    A(f"  seas in the storm {4 * np.sqrt(_m0(r.f, r.S_source)):.1f} m")
    A("")

    A(f"CROSSING  {s.distance_km:.0f} km")
    # "Fills in" at a quarter of the peak height. A 5% threshold is not useful:
    # the arrival window is Gaussian, so its tail reaches back further than any
    # real swell would, and for a nearby storm it clips against t = 0.
    A(f"  swell fills in    {_hours(r.t_hours[np.argmax(r.Hm0_t > 0.25 * r.Hm0_t.max())])}")
    A(f"  peak arrives      {_hours(p['time_h'])}")
    A(f"  chirp             {p['chirp_hz_per_day']:.4f} Hz/day"
      f"  ({p['chirp_hz_per_day'] * 1e3:.1f} mHz/day)")
    A("")

    A("AT THE BUOY")
    A(f"  Hm0               {p['Hm0']:.2f} m")
    A(f"  peak period       {p['Tp']:.1f} s")
    A(f"  bandwidth nu      {p['nu']:.3f}"
      f"   ({'clean swell' if p['nu'] < 0.15 else 'mixed' if p['nu'] < 0.3 else 'windy'})")
    A(f"  biggest in an hour ~{r.record['H_hundredth']:.2f} m")
    A(f"  BFI               {r.surf['BFI']:.2f}"
      f"   ({'Rayleigh holds' if r.surf['BFI'] < 1 else 'expect a fat tail'})")
    A("")

    A("SETS")
    A(f"  envelope corr.    kappa = {g['kappa']:.3f}")
    A(f"  waves per set     {g['mean_run_length']:.1f}  (predicted)")
    if m.get("n_waves"):
        A(f"                    {m['mean_run_length']:.1f}  (counted in the record,"
          f" {m['n_groups']} sets in {m['n_waves']} waves)")
    A(f"  set interval      {g['set_interval'] / 60:.1f} min")
    A("")

    A("SURF")
    if b is None:
        A("  the wave never breaks on this profile")
    else:
        A(f"  breaks in         {b.h:.1f} m of water, {b.x:.0f} m out")
        A(f"  breaker height    {b.H:.2f} m   (gamma_b = {b.gamma_b:.2f})")
        A(f"  Iribarren xi_0    {b.xi0:.2f}  ->  {b.kind.upper()}")
        A(f"  peel angle        {b.peel_deg:.0f} deg"
          f"   ({'closeout' if b.peel_deg < 3 else 'peels'})")
    return "\n".join(L)


def _m0(f, S):
    from .util import trapz
    return float(trapz(S, f))


def _hours(h):
    h = float(h)
    d, rem = divmod(h, 24.0)
    return f"{h:.1f} h ({d:.0f} d {rem:.0f} h)"
