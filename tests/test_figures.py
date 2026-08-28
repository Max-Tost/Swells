"""The figure script must at least import and expose every figure it claims.

It sits outside the package and manipulates sys.path to find `swells`, which
makes it easy to break with an import placed above that line -- exactly what
happened once. Nothing else in the suite touches it, so the lessons would have
gone on referencing PNGs that could no longer be regenerated.
"""

import importlib.util
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "figures", "make_figures.py")


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("make_figures", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)          # fails loudly on a bad import order
    return m


def test_every_figure_referenced_by_a_lesson_exists(mod):
    """Follow the image links in the lessons and check the files are there."""
    import glob
    import re

    referenced = set()
    for path in glob.glob(os.path.join(ROOT, "lessons", "*.md")):
        with open(path) as fh:
            referenced |= set(re.findall(r"figures/(fig[\w.]+\.png)", fh.read()))

    assert referenced, "no figures referenced by any lesson"
    missing = [f for f in sorted(referenced)
               if not os.path.exists(os.path.join(ROOT, "figures", f))]
    assert not missing, f"lessons reference missing figures: {missing}"


def test_script_defines_a_figure_function_for_each_png(mod):
    fns = [n for n in dir(mod) if n.startswith("fig_")]
    assert len(fns) >= 10
