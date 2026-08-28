"""The whole chain, end to end.

Walked through in lessons/11-the-whole-chain.md.

Storm -> spectrum -> ocean crossing -> arrival -> groups -> surf zone.
Set the dials, get everything.
"""

from dataclasses import dataclass, field

import numpy as np

from . import groups, nearshore, nonlinear, propagate, spectra, synth
from .constants import G
from .dispersion import deep_wavelength
from .fetch import sea_state


@dataclass
class Storm:
    U10: float = 25.0             # wind speed, m/s
    fetch_km: float = 600.0       # fetch, km
    duration_h: float = 36.0      # how long it blows, hours
    distance_km: float = 4000.0   # great-circle distance to the coast, km
    gamma: float = 3.3

    @property
    def fetch(self):
        return self.fetch_km * 1e3

    @property
    def duration(self):
        return self.duration_h * 3600.0

    @property
    def distance(self):
        return self.distance_km * 1e3


@dataclass
class Coast:
    slope: float = 0.02           # beach slope tan(beta)
    approach_deg: float = 30.0    # deep-water angle between crests and contours
    x_offshore: float = 2000.0


@dataclass
class Result:
    storm: Storm
    coast: Coast
    sea: object
    f: np.ndarray
    S_source: np.ndarray
    t_hours: np.ndarray
    S_ft: np.ndarray
    Hm0_t: np.ndarray
    Tp_t: np.ndarray
    peak: dict = field(default_factory=dict)
    record: dict = field(default_factory=dict)
    grouping: dict = field(default_factory=dict)
    measured: dict = field(default_factory=dict)
    surf: dict = field(default_factory=dict)
    breaking: object = None


def run(storm=None, coast=None, f=None, n_times=400, record_hours=1.0,
        dt=0.5, seed=0):
    """Run the full pipeline and return everything."""
    storm = storm or Storm()
    coast = coast or Coast()

    if f is None:
        f = np.linspace(0.015, 0.30, 1200)

    # 1. Generation ---------------------------------------------------------
    S_source, sea = spectra.storm_spectrum(
        f, storm.U10, storm.fetch, storm.duration, gamma=storm.gamma)

    # 2. Crossing -----------------------------------------------------------
    # Open the window well before the fastest component could possibly arrive,
    # so the event is fully contained and "first arrivals" means what it says.
    lead = 3.0 * propagate.arrival_duration(f.min(), storm.duration, storm.fetch)
    t_first = propagate.travel_time(f.min(), storm.distance)
    t_last = (propagate.travel_time(f.max(), storm.distance)
              + 2.0 * propagate.arrival_duration(f.max(), storm.duration, storm.fetch))
    t = np.linspace(max(t_first - lead, 0.0), t_last, n_times)
    S_ft = propagate.spectrogram(f, S_source, storm.distance,
                                 storm.duration, storm.fetch, t)
    Hm0_t, Tp_t = propagate.sea_state_at(f, S_ft)

    # 3. Peak of the event --------------------------------------------------
    i_peak = int(np.argmax(Hm0_t))
    S_peak = S_ft[i_peak]
    Hm0_peak = float(Hm0_t[i_peak])
    Tp_peak = float(Tp_t[i_peak])

    peak = {
        "time_h": float(t[i_peak] / 3600.0),
        "Hm0": Hm0_peak,
        "Tp": Tp_peak,
        "nu": spectra.bandwidth_nu(f, S_peak),
        "T_mean": spectra.mean_period(f, S_peak),
        "chirp_hz_per_day": propagate.chirp_rate(storm.distance) * 86400.0,
    }

    # 4. A real sea surface, so we can watch the sets -----------------------
    grouping = groups.group_statistics(f, S_peak)
    t_rec, eta = synth.surface_elevation(f, S_peak, record_hours * 3600.0, dt,
                                         seed=seed)
    env = groups.envelope(eta)
    measured = groups.measured_group_statistics(t_rec, eta)
    record = {"t": t_rec, "eta": eta, "envelope": env,
              **groups.height_statistics(Hm0_peak)}

    # 5. Is it nonlinear? ---------------------------------------------------
    k_p = (2.0 * np.pi * (1.0 / Tp_peak)) ** 2 / G if Tp_peak > 0 else np.nan
    m0 = spectra.moment(f, S_peak, 0)
    bfi = nonlinear.benjamin_feir_index(k_p, m0, peak["nu"])

    # 6. The surf zone ------------------------------------------------------
    x, h = nearshore.plane_beach(coast.slope, coast.x_offshore)
    table, event = nearshore.transform_transect(
        Tp_peak, Hm0_peak, coast.approach_deg, x, h, tan_beta=coast.slope)

    surf = {"table": table, "L0": deep_wavelength(1.0 / Tp_peak), "BFI": bfi}

    return Result(storm=storm, coast=coast, sea=sea, f=f, S_source=S_source,
                  t_hours=t / 3600.0, S_ft=S_ft, Hm0_t=Hm0_t, Tp_t=Tp_t,
                  peak=peak, record=record, grouping=grouping,
                  measured=measured, surf=surf, breaking=event)
