"""Cross-check with ESA's independently computed orbit (NEO Coordination Centre).

1. Carry the JPL orbit to ESA's epoch with this model and compare element by element.
2. Integrate ESA's orbit and summarise its resonant behaviour, to see whether the result
   depends on whose orbit solution is used.

    python scripts/compare_esa.py 2026_LO4
"""
import _bootstrap
import argparse

from coorb import classify, nbody, orbits, summary
from coorb.paths import object_dir, year_of_jd

NAMES = [("a_au", "a"), ("e", "e"), ("i_deg", "i"), ("node_deg", "node"), ("peri_deg", "peri"), ("M_deg", "M")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("object")
    ap.add_argument("--years", type=float, help="span each way for the ESA-orbit run")
    ap.add_argument("--dt", type=float, help="output cadence in years")
    args = ap.parse_args()
    obj = orbits.load_object(args.object)
    jpl, esa, an = obj["jpl"], obj["esa"], obj["analysis"]
    ep_j, ep_e = float(jpl["epoch_jd_tdb"]), float(esa["epoch_jd"])
    model = nbody.Model(ep_j, [orbits.element_tuple(jpl["elements"])], gr=True, split_moon=True)
    model.advance(ep_e - ep_j)
    mine = model.helio_elements(0)
    rows = []
    print(f"{obj['designation']}: JPL orbit carried {ep_e - ep_j:+.0f} days to ESA's epoch (JD {ep_e}), against ESA's orbit")
    for key, short in NAMES:
        d = mine[short] - esa["elements"][key]
        if short in ("node", "peri", "M"):
            d = (d + 180) % 360 - 180
        rows.append(dict(element=short, jpl_propagated=mine[short], esa=esa["elements"][key], diff=d, diff_in_esa_sigma=d / esa["sigma"][key]))
        print(f"   {short:<5s} JPL {mine[short]:.10f}   ESA {esa['elements'][key]:.10f}   diff {d:+.2e}  ({d / esa['sigma'][key]:+.2f} ESA sigma)")
    years = args.years or an["esa_years"]
    dt = args.dt or an["nominal_dt"]
    out = nbody.integrate_both_parallel(ep_e, [orbits.element_tuple(esa["elements"])], obj["planet"], years, dt_out_years=dt, gr=True)
    n = classify.samples_for(an["smooth_years"], dt)
    s = summary.summarise(out["t"], out["phi"][:, 0], out["a"][:, 0], out["a_pl"], out["d"][:, 0], year_of_jd(ep_e), n)
    lane = classify.lane_switches(out["t"], out["a"][:, 0], out["a_pl"], classify.samples_for(10.0, dt))
    s["lane_switch_years"] = [float(year_of_jd(ep_e) + x) for x in lane]
    lib, hs = s["librating"], s["pure_horseshoe"]
    print(f"ESA orbit, +-{years:.0f} yr: librating from {lib['start']:.0f} to {lib['end']:.0f} yr (period {lib['period']:.0f} yr); "
          f"does not cross the planet's position between {hs['start']:.0f} and {hs['end']:.0f} yr")
    print("   turning points near now:", ", ".join(f"{p['year']:.0f} ({p['deg_ahead_of_planet']:.1f} deg)" for p in s["turning_points_near_now"]))
    if "beside_planet" in s:
        print(f"   beside the planet (|phi| < {s['beside_planet']['within_deg']:.0f} deg) from {s['beside_planet']['from_year']:.0f} to {s['beside_planet']['to_year']:.0f}")
    _bootstrap.write_json(object_dir(args.object) / "esa_comparison.json",
                          dict(object=obj["designation"], esa_epoch_jd=ep_e, elements=rows, esa_orbit_run_years=years, esa_orbit=s))


if __name__ == "__main__":
    main()
