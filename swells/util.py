"""Small helpers."""

import numpy as np

# numpy 2.0 renamed trapz -> trapezoid. Support both.
trapz = getattr(np, "trapezoid", None) or np.trapz


def as_array(x):
    return np.asarray(x, dtype=np.float64)
