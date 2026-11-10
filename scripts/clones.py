"""How much does the conclusion depend on the orbit's uncertainty?

Draws clone orbits from the JPL covariance matrix, integrates the nominal orbit and every
clone (Sun, eight planets, relativity), and records for each one how long the libration and
the pure horseshoe last. Backward and forward runs use one process each.

    python scripts/clones.py 2026_LO4                         # as in the write-up (slow)
    python scripts/clones.py 2026_LO4 --n 8 --years 2000      # quick check
    python scripts/clones.py 2026_LO4 --summarise-only        # redo the tables from runs/
"""
import _bootstrap
import argparse
import csv
import time

import numpy as np

from coorb import classify, nbody, orbits, summary
from coorb.paths import RUNS, object_dir, year_of_jd


def summarise(name, obj, tag, smooth_years):
    d = np.load(RUNS / f"{name}_clones_{tag}.npz")
    t, phi = d["t"], d["phi"]
    dt = float(t[1] - t[0])
    n = classify.samples_for(smooth_years, dt)
    rows = []
    for j in range(phi.shape[1]):
        lib = classify.bounded_interval(t, phi[:, j], n)
        hs = classify.pure_horseshoe_interval(t, phi[:, j], n)
        rows.append(dict(clone=j, a_au=d["el"][j][0], e=d["el"][j][1],
                         librating_from_yr=lib["start"], librating_to_yr=lib["end"], libration_period_yr=lib["period"],
                         librating_at_run_start=lib["at_run_start"], librating_at_run_end=lib["at_run_end"],
                         horseshoe_from_yr=hs["start"], horseshoe_to_yr=hs["end"],
                         horseshoe_at_run_start=hs["at_run_start"], horseshoe_at_run_end=hs["at_run_end"]))
    path = object_dir(name) / f"clones_{tag}.csv"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.10f}" if k in ("a_au", "e") else f"{v:.1f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {_bootstrap.show(path)}")
    col = lambda k: np.array([r[k] for r in rows], dtype=float)

    def stats(x):
        ok = x[np.isfinite(x)]            # a quantity can be undefined for every orbit of a short run
        if not len(ok):
            return dict(nominal=float(x[0]), min=float("nan"), median=float("nan"), max=float("nan"))
        return dict(nominal=float(x[0]), min=float(ok.min()), median=float(np.median(ok)), max=float(ok.max()))

    i0 = int(np.argmin(np.abs(t)))
    spread = {}
    for T in (100, 500, 1000, 2000, 5000, 10000, 20000):
        if T <= t[-1]:
            spread[f"-{T}"] = summary.circular_spread_deg(phi[int(np.argmin(np.abs(t + T)))])
            spread[f"+{T}"] = summary.circular_spread_deg(phi[int(np.argmin(np.abs(t - T)))])
    s = dict(object=obj["designation"], n_orbits=len(rows), run_from_yr=float(t[0]), run_to_yr=float(t[-1]),
             epoch_year=float(year_of_jd(float(obj["jpl"]["covariance"]["epoch_jd_tdb"]))),
             phi_now_deg=[float(phi[i0].min()), float(phi[i0].max())],
             librating_from_yr=stats(col("librating_from_yr")), librating_to_yr=stats(col("librating_to_yr")),
             libration_period_yr=stats(col("libration_period_yr")),
             n_librating_at_run_start=int(col("librating_at_run_start").sum()), n_librating_at_run_end=int(col("librating_at_run_end").sum()),
             horseshoe_from_yr=stats(col("horseshoe_from_yr")), horseshoe_to_yr=stats(col("horseshoe_to_yr")),
             n_horseshoe_at_run_start=int(col("horseshoe_at_run_start").sum()), n_horseshoe_at_run_end=int(col("horseshoe_at_run_end").sum()),
             n_horseshoe_whole_run=int((col("horseshoe_at_run_start") * col("horseshoe_at_run_end")).sum()),
             spread_of_phi_among_orbits_deg=spread,
             min_sampled_distance_to_planet_au=float(np.nanmin(d["d"])))
    print(f"{obj['designation']}: {s['n_orbits']} orbits (nominal + clones), {t[0]:.0f} to {t[-1]:.0f} yr")
    for key, label in (("librating", "librating"), ("horseshoe", "pure horseshoe")):
        a, b = s[f"{key}_from_yr"], s[f"{key}_to_yr"]
        print(f"   {label:<15s} starts: nominal {a['nominal']:.0f}, range {a['min']:.0f} .. {a['max']:.0f} yr | "
              f"ends: nominal {b['nominal']:.0f}, range {b['min']:.0f} .. {b['max']:.0f} yr")
    print(f"   still librating at run start/end: {s['n_librating_at_run_start']}/{s['n_librating_at_run_end']} of {s['n_orbits']}; "
          f"pure horseshoe for the whole run: {s['n_horseshoe_whole_run']}")
    print("   spread of phi among orbits (deg):", ", ".join(f"{k} yr: {v:.1f}" for k, v in spread.items()))
    _bootstrap.write_json(object_dir(name) / f"clones_{tag}_summary.json", s)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("object")
    ap.add_argument("--n", type=int, help="number of clones (the nominal orbit is always added)")
    ap.add_argument("--years", type=float, help="span each way")
    ap.add_argument("--dt", type=float, help="output cadence in years")
    ap.add_argument("--seed", type=int)
    ap.add_argument("--tag", help="label for the output files (default: n<clones>_<years>yr)")
    ap.add_argument("--summarise-only", action="store_true", help="skip the integration; reuse runs/<object>_clones_<tag>.npz")
    args = ap.parse_args()
    obj = orbits.load_object(args.object)
    an = obj["analysis"]
    n = args.n or an["clone_n"]
    years = args.years or an["clone_years"]
    dt = args.dt or an["clone_dt"]
    seed = an["clone_seed"] if args.seed is None else args.seed
    tag = args.tag or f"n{n}_{years:.0f}yr"
    if not args.summarise_only:
        epoch, asts = orbits.sample_clones(obj, n, seed)
        t0 = time.time()
        out = nbody.integrate_both_parallel(epoch, asts, obj["planet"], years, dt_out_years=dt, gr=True)
        print(f"integrated {len(asts)} orbits, +-{years:.0f} yr, in {time.time() - t0:.0f} s")
        RUNS.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(RUNS / f"{args.object}_clones_{tag}.npz", t=out["t"], phi=out["phi"].astype(np.float32),
                            a=out["a"].astype(np.float32), e=out["e"].astype(np.float32), inc=out["inc"].astype(np.float32),
                            d=out["d"].astype(np.float32), a_pl=out["a_pl"], el=np.array(asts))
    summarise(args.object, obj, tag, an["smooth_years"])


if __name__ == "__main__":
    main()
