"""The browser app's adapter must stay JSON-serialisable and complete.

The app fetches swells/*.py at runtime, so a rename in the package silently
breaks the UI and nothing else notices. This test is the tripwire.
"""

import json
import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "app"))

import bridge  # noqa: E402


@pytest.fixture(scope="module")
def payload():
    return json.loads(bridge.simulate(25, 600, 36, 4000, 0.02, 30,
                                      record_hours=0.4))


def test_all_sections_present(payload):
    for key in ["report", "storm", "spectrum", "spectrogram", "timeline",
                "record", "buoy", "surf"]:
        assert key in payload


def test_everything_is_finite(payload):
    """A NaN reaching the canvas produces a blank plot and no error message."""
    def walk(node, path=""):
        if isinstance(node, dict):
            for k, v in node.items():
                walk(v, f"{path}.{k}")
        elif isinstance(node, list):
            for i, v in enumerate(node[:5000]):
                walk(v, f"{path}[{i}]")
        elif isinstance(node, float):
            assert math.isfinite(node), f"non-finite at {path}"
    walk(payload)


def test_array_lengths_match(payload):
    s = payload["spectrum"]
    assert len(s["f"]) == len(s["S_source"]) == len(s["S_peak"])
    r = payload["record"]
    assert len(r["t"]) == len(r["eta"]) == len(r["env"])
    sg = payload["spectrogram"]
    assert len(sg["z"]) == len(sg["t"])
    assert all(len(row) == len(sg["f"]) for row in sg["z"])
    su = payload["surf"]
    assert len(su["x"]) == len(su["h"]) == len(su["H"])


def test_downsampling_stays_within_canvas_budget(payload):
    assert len(payload["record"]["t"]) <= 2400
    assert len(payload["spectrogram"]["t"]) <= 220
    assert len(payload["spectrogram"]["f"]) <= 140


def test_spectrogram_has_energy_in_it(payload):
    """Cropping the axes must not crop the swell out of the picture."""
    assert payload["spectrogram"]["zmax"] > 0
    z = payload["spectrogram"]["z"]
    peak_row = max(range(len(z)), key=lambda i: max(z[i]))
    assert 0 < peak_row < len(z) - 1      # the event is inside the window


def test_extremes_do_not_crash():
    """Slider endpoints, in the combinations most likely to break something."""
    for args in [(8, 50, 3, 200, 1 / 200, 0),      # everything minimal
                 (40, 2500, 96, 16000, 1 / 5, 70),  # everything maximal
                 (40, 50, 96, 200, 1 / 5, 0),     # short fetch, close, steep
                 (8, 2500, 3, 16000, 1 / 200, 70)]:  # weak, brief, far
        d = json.loads(bridge.simulate(*args, record_hours=0.25))
        assert d["buoy"]["Hm0"] >= 0
        assert d["storm"]["Tp"] > 0


def test_no_numpy2_only_apis_in_pyodide_path():
    """Pyodide ships numpy 1.26, which has np.trapz but not np.trapezoid.

    Anything the browser imports -- the whole swells package plus bridge.py --
    must therefore go through swells.util.trapz rather than calling the numpy 2
    spelling directly. This exact bug shipped once and only showed up as a
    traceback in the browser console, which no local test would have caught.
    """
    import glob
    import re

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    banned = re.compile(r"\bnp\.trapezoid\b")
    offenders = []
    for path in glob.glob(os.path.join(root, "swells", "*.py")) + \
            [os.path.join(root, "app", "bridge.py")]:
        if os.path.basename(path) == "util.py":
            continue          # the shim itself is allowed to know both names
        with open(path) as fh:
            for n, line in enumerate(fh, 1):
                if banned.search(line):
                    offenders.append(f"{os.path.relpath(path, root)}:{n}")
    assert not offenders, ("numpy-2-only API on the Pyodide import path: "
                           + ", ".join(offenders))
