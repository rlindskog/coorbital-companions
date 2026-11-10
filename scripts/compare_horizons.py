"""Compare this model with NASA JPL Horizons' own integration of the same orbit.

data/horizons/<object>.csv holds osculating elements from Horizons every 25 years from 1600
to 2500. This script integrates the same JPL orbit solution with the model here (Sun, eight
planets, Moon, relativity) and compares the two at every date.

    python scripts/compare_horizons.py 2026_LO4
"""
import _bootstrap
import argparse
import csv

import numpy as np

from coorb import nbody, orbits
from coorb.ephem import GM
from coorb.paths import object_dir, year_of_jd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("object")
    args = ap.parse_args()
    obj = orbits.load_object(args.object)
    H = orbits.load_horizons(args.object)
    ep0 = float(obj["jpl"]["epoch_jd_tdb"])
    el = orbits.element_tuple(obj["jpl"]["elements"])
    planet = obj["planet"]
    mu = GM["sun"]
    res = {}
    for direction in (-1, +1):
        model = nbody.Model(ep0, [el], gr=True, split_moon=True)
        rows = H[H[:, 0] < ep0][::-1] if direction < 0 else H[H[:, 0] > ep0]
        for row in rows:
            buf = model.advance(row[0] - ep0)
            sun, ast, pl = buf[0], buf[model.n_massive], model.planet_state(planet)
            a, e, inc, lam = nbody.elements(mu, ast[:3] - sun[:3], ast[3:] - sun[3:])
            ap_, ep_, ip_, lp = nbody.elements(mu, pl[:3] - sun[:3], pl[3:] - sun[3:])
            res[row[0]] = (float(a), float(e), float(np.degrees(lam)), float(np.degrees(lp)), float(ap_))
    table = []
    for jd, ec, inc, node, peri, ma, a in H:
        lam_h = (node + peri + ma) % 360
        a_m, e_m, lam_m, lam_p, a_p = res[jd]
        table.append(dict(year=round(year_of_jd(jd), 3), jd=jd, a_horizons_au=a, a_model_au=a_m, e_horizons=ec, e_model=e_m,
                          dlambda_deg=(lam_m - lam_h + 180) % 360 - 180, deg_ahead_of_planet=(lam_h - lam_p) % 360,
                          a_planet_au=a_p))
    out = object_dir(args.object) / "horizons_comparison.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(table[0]))
        w.writeheader()
        for r in table:
            w.writerow({k: (f"{v:.10f}" if k.startswith(("a_", "e_")) else f"{v:.6f}" if k in ("dlambda_deg", "deg_ahead_of_planet") else v) for k, v in r.items()})
    print(f"wrote {_bootstrap.show(out)}")
    da = max(abs(r["a_model_au"] - r["a_horizons_au"]) for r in table)
    dl = max(abs(r["dlambda_deg"]) for r in table)
    de = max(abs(r["e_model"] - r["e_horizons"]) for r in table)
    print(f"{obj['designation']} vs JPL Horizons ({obj['horizons']['solution']}), {len(table)} dates {table[0]['year']:.0f}-{table[-1]['year']:.0f}:")
    print(f"   largest difference in mean longitude {dl:.6f} deg, semi-major axis {da:.1e} au, eccentricity {de:.1e}")
    print("   year   deg ahead of planet   a - a_planet (au)   [from the Horizons elements]")
    for r in table[::4]:
        print(f"   {r['year']:6.0f}   {r['deg_ahead_of_planet']:8.2f}            {r['a_horizons_au'] - r['a_planet_au']:+.5f}")
    _bootstrap.write_json(object_dir(args.object) / "horizons_summary.json",
                          dict(object=obj["designation"], horizons_solution=obj["horizons"]["solution"], n_dates=len(table),
                               first_year=table[0]["year"], last_year=table[-1]["year"], max_abs_dlambda_deg=dl,
                               max_abs_da_au=da, max_abs_de=de))


if __name__ == "__main__":
    main()
