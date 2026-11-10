"""Turn one simulated time series into the handful of numbers quoted in the write-up."""
import numpy as np

from . import classify


def summarise(t, phi, a, a_pl, d, epoch_year, smooth_n, prom=30.0, near_deg=25.0):
    """Summary of one asteroid's resonant behaviour. `t` is in years from the epoch."""
    t = np.asarray(t, dtype=float)
    phi = np.asarray(phi, dtype=float)
    i0 = int(np.argmin(np.abs(t)))
    tp, pp, ps = classify.turning_points(t, phi, smooth_n, prom)
    out = dict(epoch_year=float(epoch_year), run_from_yr=float(t[0]), run_to_yr=float(t[-1]),
               phi_now_deg=float(phi[i0]), deg_ahead_of_planet_now=float(np.mod(phi[i0], 360.0)))
    out["librating"] = classify.bounded_interval(t, phi, smooth_n, prom)
    out["pure_horseshoe"] = classify.pure_horseshoe_interval(t, phi, smooth_n)
    order = np.argsort(np.abs(tp))[:6]
    out["turning_points_near_now"] = [dict(year=float(epoch_year + tp[j]), deg_ahead_of_planet=float(np.mod(pp[j], 360.0)))
                                      for j in sorted(order, key=lambda j: tp[j])]
    da = np.asarray(a) - np.asarray(a_pl)
    out["a_minus_a_planet_au"] = [float(np.nanmin(da)), float(np.nanmax(da))]
    k = int(np.nanargmin(d))
    out["min_sampled_distance_to_planet"] = dict(au=float(d[k]), year=float(epoch_year + t[k]))
    if abs(phi[i0]) < near_deg:                       # currently beside the planet: how long?
        near = np.abs(phi) < near_deg
        lo = i0
        while lo > 0 and near[lo - 1]:
            lo -= 1
        hi = i0
        while hi < len(t) - 1 and near[hi + 1]:
            hi += 1
        out["beside_planet"] = dict(within_deg=near_deg, from_year=float(epoch_year + t[lo]), to_year=float(epoch_year + t[hi]))
    return out


def circular_spread_deg(phi_row):
    """Circular standard deviation (degrees) of a set of angles."""
    z = np.exp(1j * np.deg2rad(np.asarray(phi_row, dtype=float)))
    r = min(1.0, abs(z.mean()))
    return float(np.degrees(np.sqrt(-2 * np.log(r)))) if r > 0 else float("inf")
