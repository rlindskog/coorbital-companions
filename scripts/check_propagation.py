"""Check the force model against JPL's own propagation of the orbit.

JPL publishes the orbit twice: osculating elements at a standard epoch, and (with the
covariance) at the orbit-solution epoch. Integrating one to the other and comparing is an
end-to-end test of this model (frames, units, masses, relativity) and of the input data.

    python scripts/check_propagation.py 2026_LO4
"""
import _bootstrap
import argparse

from coorb import nbody, orbits
from coorb.paths import object_dir


def run(obj, gr, split_moon):
    jpl = obj["jpl"]
    cov = jpl["covariance"]
    ep0, epc = float(jpl["epoch_jd_tdb"]), float(cov["epoch_jd_tdb"])
    model = nbody.Model(ep0, [orbits.element_tuple(jpl["elements"])], gr=gr, split_moon=split_moon)
    model.advance(epc - ep0)
    el = model.helio_elements(0)
    period = 360.0 / el["n"]
    m = el["M"] if el["M"] < 180 else el["M"] - 360
    tp = epc - m / el["n"]
    tp_ref = cov["nominal"][2]
    tp += period * round((tp_ref - tp) / period)
    mine = dict(e=el["e"], q=el["q"], tp=tp, node=el["node"], peri=el["peri"], i=el["i"])
    sig = [cov["matrix"][k][k] ** 0.5 for k in range(6)]
    rows = []
    for k, lab in enumerate(cov["labels"]):
        d = mine[lab] - cov["nominal"][k]
        rows.append(dict(element=lab, mine=mine[lab], jpl=cov["nominal"][k], diff=d, diff_sigma=d / sig[k]))
    return rows, epc - ep0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("object")
    args = ap.parse_args()
    obj = orbits.load_object(args.object)
    out = {"object": obj["designation"]}
    for label, gr, split in [("newtonian", False, True), ("with_relativity", True, True)]:
        rows, span = run(obj, gr, split)
        out[label] = rows
        out["span_days"] = span
        print(f"{obj['designation']}: {label}, {abs(span):.0f} days from the element epoch to the covariance epoch")
        for r in rows:
            print(f"   {r['element']:<5s} mine {r['mine']:.12f}   JPL {r['jpl']:.12f}   diff {r['diff']:+.2e}  ({r['diff_sigma']:+.4f} sigma)")
    _bootstrap.write_json(object_dir(args.object) / "propagation_check.json", out)


if __name__ == "__main__":
    main()
