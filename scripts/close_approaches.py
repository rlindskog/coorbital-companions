"""Actual distances from the asteroid to its planet and to Earth, at fine time resolution.

The other scripts sample a few times per year, which is too coarse to catch the closest
passes. This one steps every day or two over a limited span (set in the object file).

    python scripts/close_approaches.py 2026_LO4
"""
import _bootstrap
import argparse

import numpy as np

from coorb import nbody, orbits
from coorb.paths import object_dir, year_of_jd

AU_MKM = 149.5978707      # million km per au


def one_direction(obj, direction, years, step):
    ep0 = float(obj["jpl"]["epoch_jd_tdb"])
    model = nbody.Model(ep0, [orbits.element_tuple(obj["jpl"]["elements"])], gr=True, split_moon=True)
    ts = direction * np.arange(0, years * 365.25 + step, step)
    out = np.zeros((len(ts), 3))
    ie = model.index("earth")
    for k, t in enumerate(ts):
        buf = model.advance(t)
        ast = buf[model.n_massive, :3]
        out[k] = [t, np.linalg.norm(ast - model.planet_state(obj["planet"])[:3]), np.linalg.norm(ast - buf[ie, :3])]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("object")
    ap.add_argument("--from-year", type=float)
    ap.add_argument("--to-year", type=float)
    ap.add_argument("--step-days", type=float)
    ap.add_argument("--list-below-au", type=float, default=0.12, help="list every planet approach closer than this")
    args = ap.parse_args()
    obj = orbits.load_object(args.object)
    an = obj["analysis"]
    y0 = year_of_jd(float(obj["jpl"]["epoch_jd_tdb"]))
    lo = args.from_year or an["close_from_year"]
    hi = args.to_year or an["close_to_year"]
    step = args.step_days or an["close_step_days"]
    fw = one_direction(obj, +1, hi - y0, step)
    bw = one_direction(obj, -1, y0 - lo, step)
    allr = np.vstack([bw[::-1][:-1], fw])
    yr = y0 + allr[:, 0] / 365.25
    dp, de = allr[:, 1], allr[:, 2]
    i, j = int(np.argmin(dp)), int(np.argmin(de))
    planet = obj["planet"].capitalize()
    print(f"{obj['designation']}, {lo:.0f}-{hi:.0f}, sampled every {step:g} d")
    print(f"   closest to {planet}: {dp[i]:.4f} au ({dp[i] * AU_MKM:.1f} million km) in {yr[i]:.1f}")
    print(f"   closest to Earth: {de[j]:.4f} au ({de[j] * AU_MKM:.1f} million km) in {yr[j]:.1f}")
    loc = [k for k in range(1, len(dp) - 1) if dp[k] < dp[k - 1] and dp[k] <= dp[k + 1] and dp[k] < args.list_below_au]
    groups = []                                   # cluster approaches separated by less than 5 years
    for k in loc:
        if groups and yr[k] - yr[groups[-1][-1]] < 5:
            groups[-1].append(k)
        else:
            groups.append([k])
    clusters = []
    for g in groups:
        kmin = min(g, key=lambda k: dp[k])
        clusters.append(dict(from_year=float(yr[g[0]]), to_year=float(yr[g[-1]]), n_passes=len(g),
                             closest_au=float(dp[kmin]), closest_year=float(yr[kmin])))
        print(f"   {len(g):3d} passes within {args.list_below_au} au of {planet} between {yr[g[0]]:.1f} and {yr[g[-1]]:.1f}; "
              f"closest {dp[kmin]:.4f} au in {yr[kmin]:.1f}")
    _bootstrap.write_json(object_dir(args.object) / "close_approaches.json",
                          dict(object=obj["designation"], from_year=lo, to_year=hi, step_days=step,
                               closest_to_planet=dict(au=float(dp[i]), million_km=float(dp[i] * AU_MKM), year=float(yr[i])),
                               closest_to_earth=dict(au=float(de[j]), million_km=float(de[j] * AU_MKM), year=float(yr[j])),
                               approach_clusters_below_au=args.list_below_au, approach_clusters=clusters))


if __name__ == "__main__":
    main()
