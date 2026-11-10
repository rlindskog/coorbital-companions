"""Input data: object files (orbit + covariance), zone candidate lists, Horizons tables."""
import csv
import json

import numpy as np

from .ephem import GM
from .paths import DATA

ELEMENT_KEYS = ["a_au", "e", "i_deg", "node_deg", "peri_deg", "M_deg"]


def load_object(name):
    """Read data/objects/<name>.json (for example "2026_LO4")."""
    with open(DATA / "objects" / f"{name}.json") as fh:
        return json.load(fh)


def element_tuple(d):
    """(a, e, i, node, peri, M) from a dict keyed by ELEMENT_KEYS."""
    return tuple(float(d[k]) for k in ELEMENT_KEYS)


def kepler_from_cov_params(p, epoch_jd):
    """(e, q, tp, node, peri, i), the JPL covariance parameter set, -> (a, e, i, node, peri, M)."""
    e, q, tp, node, peri, inc = (float(x) for x in p)
    a = q / (1 - e)
    n = np.sqrt(GM["sun"] / a ** 3)                 # rad/day
    M = np.degrees(n * (epoch_jd - tp)) % 360
    return (a, e, inc, node, peri, M)


def sample_clones(obj, n, seed):
    """Nominal orbit followed by n clones drawn from the JPL covariance.

    Clones are multivariate-normal draws in the covariance parameters (e, q, tp, node, peri,
    i) at the covariance epoch. Returns (epoch_jd, list of element tuples).
    """
    cov = obj["jpl"]["covariance"]
    epoch = float(cov["epoch_jd_tdb"])
    nominal = np.array(cov["nominal"], dtype=float)
    C = np.array(cov["matrix"], dtype=float)
    if not np.allclose(C, C.T, rtol=1e-9, atol=0):
        raise ValueError("covariance matrix is not symmetric")
    w, v = np.linalg.eigh(C)
    L = v * np.sqrt(np.clip(w, 0, None))
    rng = np.random.default_rng(seed)
    draws = nominal + rng.standard_normal((n, 6)) @ L.T
    return epoch, [kepler_from_cov_params(nominal, epoch)] + [kepler_from_cov_params(p, epoch) for p in draws]


def load_zone(zone):
    """Rows of data/zones/<zone>.csv as dicts with numeric elements."""
    rows = []
    with open(DATA / "zones" / f"{zone}.csv") as fh:
        for r in csv.DictReader(fh):
            for k in ("epoch_jd_tdb", "a_au", "e", "i_deg", "node_deg", "peri_deg", "M_deg"):
                r[k] = float(r[k])
            rows.append(r)
    return rows


def load_horizons(name):
    """JPL Horizons osculating elements: array of rows (jd, e, i, node, peri, M, a)."""
    return np.loadtxt(DATA / "horizons" / f"{name}.csv", delimiter=",", comments="#")
