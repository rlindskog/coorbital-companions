"""Repository paths and small time helpers.

Output locations can be redirected with the environment variables COORB_RESULTS_DIR,
COORB_FIGURES_DIR and COORB_RUNS_DIR (used by `scripts/run_all.sh` in quick mode so that a
short test run does not overwrite the committed results).
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / os.environ.get("COORB_RESULTS_DIR", "results")
FIGURES = ROOT / os.environ.get("COORB_FIGURES_DIR", "figures")
RUNS = ROOT / os.environ.get("COORB_RUNS_DIR", "runs")     # large raw simulation outputs; git-ignored

J2000_JD = 2451544.5          # 2000-Jan-01 00:00


def year_of_jd(jd):
    """Decimal year (Julian years from 2000-Jan-01) of a Julian date."""
    return 2000.0 + (jd - J2000_JD) / 365.25


def object_dir(name):
    d = RESULTS / name
    d.mkdir(parents=True, exist_ok=True)
    return d
