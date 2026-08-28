"""Swells: how a storm at sea becomes sets of waves at a beach.

The lessons in lessons/ derive the physics; this package is that same physics
as code. Module by module:

    dispersion  linear wave theory, omega^2 = g k tanh(k h)      lessons 01, 02
    fetch       how big a sea the wind can build                 lesson 04
    spectra     JONSWAP, Pierson-Moskowitz, moments              lesson 04
    propagate   the ocean crossing and the dispersive sorting    lessons 05, 06
    synth       spectrum -> an actual sea surface                lesson 07
    groups      envelopes, height statistics, sets               lesson 07
    nonlinear   NLS, Benjamin-Feir, infragravity                 lesson 08
    nearshore   shoaling, refraction, set-up, breaking           lessons 09, 10
    simulate    all of it at once                                lesson 11
"""

from . import (dispersion, fetch, groups, nearshore, nonlinear, propagate,
               report, simulate, spectra, synth)
from .constants import G, RHO_A, RHO_W, R_EARTH
from .report import surf_report
from .simulate import Coast, Storm, run

__all__ = ["dispersion", "fetch", "spectra", "propagate", "synth", "groups",
           "nonlinear", "nearshore", "simulate", "report",
           "Storm", "Coast", "run", "surf_report",
           "G", "RHO_W", "RHO_A", "R_EARTH"]
