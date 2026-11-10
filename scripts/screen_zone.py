"""Screen every catalogued asteroid in a planet's co-orbital zone.

Each object in data/zones/<zone>.csv is integrated backward and forward (Newtonian, Sun +
eight planets) and its resonant angle is classified. This is the discovery step: it is quick
and approximate, and anything interesting is then examined with the per-object scripts.

    python scripts/screen_zone.py venus
    python scripts/screen_zone.py mars --years 2000        # shorter, faster
    python scripts/screen_zone.py venus --classify-only    # redo the table from runs/
"""
import _bootstrap
import argparse
import csv
import time
from multiprocessing import Pool

import numpy as np

from coorb import classify, nbody, orbits
from coorb.paths import RESULTS, RUNS

ZONES = {   # planet, default span (yr each way), output cadence (yr), smoothing window (yr)
    "venus": dict(planet="venus", years=3000, dt=0.1, smooth_years=2.0),
    "earth": dict(planet="emb", years=1500, dt=0.2, smooth_years=3.0),
    "mars": dict(planet="mars", years=6000, dt=0.5, smooth_years=5.0),
}
ORDER = {"QS": 0, "L4": 1, "L5": 2, "HS": 3, "compound": 4, "edge": 5, "passing": 6, "unstable": 7}


def _work(args):
    jd, asts, planet, years, dt = args
    return nbody.integrate_both(jd, asts, planet, years, dt_out_years=dt, gr=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("zone", choices=sorted(ZONES))
    ap.add_argument("--years", type=float, help="span each way (default depends on the zone)")
    ap.add_argument("--dt", type=float, help="output cadence in years")
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--batch", type=int, default=8, help="asteroids per simulation")
    ap.add_argument("--classify-only", action="store_true", help="skip the integration; reuse runs/screening_<zone>.npz")
    args = ap.parse_args()
    cfg = dict(ZONES[args.zone])
    rows = orbits.load_zone(args.zone)
    series_file = RUNS / f"screening_{args.zone}.npz"
    if args.classify_only:
        store = dict(np.load(series_file))
    else:
        years = args.years or cfg["years"]
        dt = args.dt or cfg["dt"]
        jobs, index = [], []
        for ep in sorted({r["epoch_jd_tdb"] for r in rows}):
            grp = [r for r in rows if r["epoch_jd_tdb"] == ep]
            for k in range(0, len(grp), args.batch):
                b = grp[k:k + args.batch]
                jobs.append((ep, [orbits.element_tuple(r) for r in b], cfg["planet"], years, dt))
                index.append(b)
        t0 = time.time()
        with Pool(args.procs) as pool:
            outs = pool.map(_work, jobs, chunksize=1)
        print(f"integrated {len(rows)} objects, +-{years:.0f} yr, in {time.time() - t0:.0f} s")
        store = {}                      # per object: rows t [yr], phi [deg], a [au], a_planet [au]
        for b, out in zip(index, outs):
            for j, r in enumerate(b):
                store[r["designation"]] = np.vstack([out["t"], out["phi"][:, j], out["a"][:, j], out["a_pl"]]).astype(np.float32)
        RUNS.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(series_file, **store)
    table = []
    for r in rows:
        t, phi = store[r["designation"]][0].astype(float), store[r["designation"]][1].astype(float)
        years, dt = float(t[-1]), float(t[1] - t[0])
        c = classify.classify_now(t, phi, smooth_n=classify.samples_for(cfg["smooth_years"], dt), prom=8.0)
        table.append(dict(designation=r["designation"], state=c["state"], phi_now_deg=c.get("phi0"),
                          turn_lo_deg=c.get("lo"), turn_hi_deg=c.get("hi"), librating_from_yr=c.get("t_start"),
                          librating_to_yr=c.get("t_end"), span_yr=c.get("span", 0.0),
                          condition_code=r["condition_code"], data_arc_days=r["data_arc_days"],
                          first_obs=r["first_obs"], a_au=r["a_au"], e=r["e"], i_deg=r["i_deg"]))
    table.sort(key=lambda c: (ORDER.get(c["state"], 9), c["designation"]))
    out_csv = RESULTS / "screening" / f"{args.zone}.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0]))
        w.writeheader()
        for c in table:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) and k not in ("a_au", "e", "i_deg") else v) for k, v in c.items()})
    print(f"wrote {_bootstrap.show(out_csv)}  (span +-{years:.0f} yr)")
    counts = {}
    for c in table:
        counts[c["state"]] = counts.get(c["state"], 0) + 1
    print("states now:", ", ".join(f"{k} {v}" for k, v in sorted(counts.items(), key=lambda kv: ORDER.get(kv[0], 9))))
    for c in table:
        if c["state"] in ("passing", "unstable"):
            continue
        fmt = lambda v: "   n/a" if v is None else f"{v:7.1f}"
        print(f"  {c['designation']:<11s} {c['state']:<8s} phi now {fmt(c['phi_now_deg'])}  turns at {fmt(c['turn_lo_deg'])} /{fmt(c['turn_hi_deg'])}"
              f"  librating {fmt(c['librating_from_yr'])} ..{fmt(c['librating_to_yr'])} yr  U={c['condition_code']} arc={c['data_arc_days'] or '?'} d")


if __name__ == "__main__":
    main()
