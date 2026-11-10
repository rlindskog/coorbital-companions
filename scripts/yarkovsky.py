"""Would a slow thermal drift (the Yarkovsky effect) break the resonance?

Small asteroids drift in semi-major axis because they re-radiate sunlight unevenly. The
drift rate of these objects is unknown, so this script brackets it: the nominal orbit is
integrated with several constant da/dt values (listed in the object file; ~1e-3 au/Myr is
typical for a 100 m body near Mars) and the lifetime of the horseshoe is recorded for each.

    python scripts/yarkovsky.py 2025_ON9
"""
import _bootstrap
import argparse
import csv
import time

import numpy as np

from coorb import classify, nbody, orbits
from coorb.paths import RUNS, object_dir


def summarise(name, obj):
    d = np.load(RUNS / f"{name}_yarkovsky.npz")
    t, phi, rates = d["t"], d["phi"], d["rates"]
    n = classify.samples_for(obj["analysis"]["smooth_years"], float(t[1] - t[0]))
    path = object_dir(name) / "yarkovsky.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["drift_au_per_myr", "librating_from_yr", "librating_to_yr", "horseshoe_from_yr", "horseshoe_to_yr",
                    "horseshoe_at_run_start", "horseshoe_at_run_end"])
        print(f"{obj['designation']}: nominal orbit with a constant drift in semi-major axis, {t[0]:.0f} to {t[-1]:.0f} yr")
        for k, r in enumerate(rates):
            lib = classify.bounded_interval(t, phi[:, k], n)
            hs = classify.pure_horseshoe_interval(t, phi[:, k], n)
            w.writerow([f"{r:+.0e}", f"{lib['start']:.0f}", f"{lib['end']:.0f}", f"{hs['start']:.0f}", f"{hs['end']:.0f}",
                        hs["at_run_start"], hs["at_run_end"]])
            print(f"   da/dt {r:+.0e} au/Myr: pure horseshoe from {hs['start']:.0f} to {hs['end']:.0f} yr")
    print(f"wrote {_bootstrap.show(path)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("object")
    ap.add_argument("--years", type=float)
    ap.add_argument("--dt", type=float, default=1.0)
    ap.add_argument("--summarise-only", action="store_true")
    args = ap.parse_args()
    obj = orbits.load_object(args.object)
    an = obj["analysis"]
    if "drift_rates_au_per_myr" not in an:
        raise SystemExit(f"no drift rates listed for {args.object} in its object file")
    if not args.summarise_only:
        rates = [float(r) for r in an["drift_rates_au_per_myr"]]
        years = args.years or an["drift_years"]
        epoch, asts = orbits.sample_clones(obj, 0, 0)          # nominal orbit at the covariance epoch
        t0 = time.time()
        out = nbody.integrate_both_parallel(epoch, [asts[0]] * len(rates), obj["planet"], years, dt_out_years=args.dt,
                                            gr=True, drift_au_per_myr=rates)
        print(f"integrated {len(rates)} drift cases, +-{years:.0f} yr, in {time.time() - t0:.0f} s")
        RUNS.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(RUNS / f"{args.object}_yarkovsky.npz", t=out["t"], phi=out["phi"], a=out["a"], a_pl=out["a_pl"],
                            rates=np.array(rates))
    summarise(args.object, obj)


if __name__ == "__main__":
    main()
