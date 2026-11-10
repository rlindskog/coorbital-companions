"""Integrate the nominal (best-fit) JPL orbit and summarise its behaviour.

Model: Sun, eight planets (Earth and Moon merged), relativity. Writes a summary and a
compact time series (used by the figure script) to results/<object>/, and the full series to
runs/.

    python scripts/nominal.py 2026_LO4
    python scripts/nominal.py 2026_LO4 --years 2000     # quick look
    python scripts/nominal.py 2026_LO4 --summarise-only # redo the summary from runs/
"""
import _bootstrap
import argparse
import csv
import time

import numpy as np

from coorb import classify, nbody, orbits, summary
from coorb.paths import RUNS, object_dir, year_of_jd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("object")
    ap.add_argument("--years", type=float, help="span each way")
    ap.add_argument("--dt", type=float, help="output cadence in years")
    ap.add_argument("--summarise-only", action="store_true", help="skip the integration; reuse runs/<object>_nominal.npz")
    args = ap.parse_args()
    obj = orbits.load_object(args.object)
    an = obj["analysis"]
    ep0 = float(obj["jpl"]["epoch_jd_tdb"])
    if args.summarise_only:
        out = dict(np.load(RUNS / f"{args.object}_nominal.npz"))
    else:
        t0 = time.time()
        out = nbody.integrate_both_parallel(ep0, [orbits.element_tuple(obj["jpl"]["elements"])], obj["planet"],
                                            args.years or an["nominal_years"], dt_out_years=args.dt or an["nominal_dt"], gr=True)
        print(f"{obj['designation']}: nominal orbit integrated in {time.time() - t0:.0f} s")
        RUNS.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(RUNS / f"{args.object}_nominal.npz", **out)
    years, dt = float(out["t"][-1]), float(out["t"][1] - out["t"][0])
    print(f"{obj['designation']}: nominal orbit, +-{years:.0f} yr")
    t, phi, a, d = out["t"], out["phi"][:, 0], out["a"][:, 0], out["d"][:, 0]
    y0 = year_of_jd(ep0)
    n = classify.samples_for(an["smooth_years"], dt)
    s = summary.summarise(t, phi, a, out["a_pl"], d, y0, n)
    s["e_range"] = [float(out["e"].min()), float(out["e"].max())]
    s["inclination_range_deg"] = [float(np.degrees(out["inc"].min())), float(np.degrees(out["inc"].max()))]
    lane = classify.lane_switches(t, a, out["a_pl"], classify.samples_for(10.0, dt))
    s["n_lane_switches"] = int(len(lane))
    lib, hs = s["librating"], s["pure_horseshoe"]
    print(f"   now {s['deg_ahead_of_planet_now']:.1f} deg ahead of the planet (phi = {s['phi_now_deg']:+.1f})")
    print(f"   librating from {lib['start']:.0f} to {lib['end']:.0f} yr relative to {y0:.1f}  (period now {lib['period']:.0f} yr)")
    print(f"   does not cross the planet's position between {hs['start']:.0f} and {hs['end']:.0f} yr "
          f"(closest {hs['gap_ahead']:.1f} deg ahead, {hs['gap_behind']:.1f} deg behind)")
    print("   turning points near now:", ", ".join(f"{p['year']:.0f} ({p['deg_ahead_of_planet']:.1f} deg)" for p in s["turning_points_near_now"]))
    if "beside_planet" in s:
        b = s["beside_planet"]
        print(f"   beside the planet (|phi| < {b['within_deg']:.0f} deg) from {b['from_year']:.0f} to {b['to_year']:.0f}")
    _bootstrap.write_json(object_dir(args.object) / "nominal_summary.json", dict(object=obj["designation"], span_years=years, **s))
    lo, hi = an.get("series_from_year", y0 - 600), an.get("series_to_year", y0 + 600)
    m = (y0 + t >= lo) & (y0 + t <= hi)
    path = object_dir(args.object) / "nominal_series.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["year", "phi_deg", "a_au", "a_planet_au", "dist_planet_au"])
        for k in np.where(m)[0]:
            w.writerow([f"{y0 + t[k]:.3f}", f"{phi[k]:.4f}", f"{a[k]:.8f}", f"{out['a_pl'][k]:.8f}", f"{d[k]:.5f}"])
    print(f"wrote {_bootstrap.show(path)}  ({int(m.sum())} rows, {lo:.0f}-{hi:.0f})")


if __name__ == "__main__":
    main()
