"""Download fresh copies of the inputs from NASA/JPL and ESA (optional).

The analysis does NOT use this script: the write-up is based on the snapshot in data/,
retrieved on 2026-11-10, and orbit solutions are revised whenever new observations arrive.
Fresh files go to data/fresh/ so the snapshot is never overwritten.

    python scripts/fetch_data.py object "2026 LO4"     # SBDB orbit + covariance, Horizons table, ESA orbit
    python scripts/fetch_data.py zone venus            # every catalogued object in the zone

CAUTION: this script was written without access to these servers (the environment in which
the repository was built could not reach them directly), so it has not been run against the
live APIs. Raw responses are always saved first; if parsing fails, start from those.
API documentation: https://ssd-api.jpl.nasa.gov/doc/sbdb.html, .../doc/sbdb_query.html,
https://ssd-api.jpl.nasa.gov/doc/horizons.html.
"""
import _bootstrap  # noqa: F401
import argparse
import csv
import json
import time

from coorb.paths import DATA

SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api"
SBDB_QUERY = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
NEOCC = "https://neo.ssa.esa.int/PSDB-portlet/download"
FRESH = DATA / "fresh"

# Semi-major-axis ranges: planet's a +- (mass ratio / 3)^(1/3) * a, the classical co-orbital width.
ZONES = {"venus": (0.7165, 0.7301), "earth": (0.990, 1.010), "mars": (1.5164, 1.5310)}


def get(url, params):
    import requests
    r = requests.get(url, params=params, timeout=120)
    r.raise_for_status()
    time.sleep(1.0)            # be polite: JPL asks for one request at a time
    return r


def fetch_zone(zone):
    lo, hi = ZONES[zone]
    cdata = json.dumps({"AND": [f"a|RG|{lo}|{hi}", "condition_code|LE|5"]})
    fields = "pdes,epoch,a,e,i,om,w,ma,condition_code,data_arc,first_obs"
    js = get(SBDB_QUERY, {"fields": fields, "sb-cdata": cdata, "sort": "a", "full-prec": "1"}).json()
    FRESH.mkdir(parents=True, exist_ok=True)
    (FRESH / f"zone_{zone}_raw.json").write_text(json.dumps(js))
    out = FRESH / f"zone_{zone}.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["designation", "epoch_jd_tdb", "a_au", "e", "i_deg", "node_deg", "peri_deg", "M_deg",
                    "condition_code", "data_arc_days", "first_obs"])
        for row in js["data"]:
            w.writerow(["" if v is None else str(v).strip() for v in row])
    print(f"{zone}: {len(js['data'])} objects -> {out}")


def fetch_object(des):
    name = des.replace(" ", "_")
    FRESH.mkdir(parents=True, exist_ok=True)
    # --- JPL Small-Body Database: orbit, covariance, absolute magnitude
    js = get(SBDB, {"sstr": des, "cov": "mat", "phys-par": "1", "full-prec": "1"}).json()
    (FRESH / f"{name}_sbdb_raw.json").write_text(json.dumps(js, indent=1))
    print(f"saved raw SBDB record -> {FRESH / (name + '_sbdb_raw.json')}")
    try:
        orb = js["orbit"]
        el = {e["name"]: float(e["value"]) for e in orb["elements"]}
        cov = orb.get("covariance") or js.get("covariance")
        cel = {e["name"]: float(e["value"]) for e in cov["elements"]}
        obj = dict(designation=des, retrieved=time.strftime("%Y-%m-%d"), jpl=dict(
            solution_id=orb.get("orbit_id"), solution_date=orb.get("soln_date"), n_obs=orb.get("n_obs_used"),
            first_obs=orb.get("first_obs"), last_obs=orb.get("last_obs"), data_arc_days=orb.get("data_arc"),
            condition_code=orb.get("condition_code"), epoch_jd_tdb=float(orb["epoch"]),
            elements=dict(a_au=el["a"], e=el["e"], i_deg=el["i"], node_deg=el["om"], peri_deg=el["w"], M_deg=el["ma"]),
            covariance=dict(epoch_jd_tdb=float(cov["epoch"]), labels=["e", "q", "tp", "node", "peri", "i"],
                            nominal=[cel["e"], cel["q"], cel["tp"], cel["om"], cel["w"], cel["i"]],
                            matrix=[[float(x) for x in row[:6]] for row in cov["data"][:6]])))
        if cov.get("labels", [])[:6] != ["e", "q", "tp", "node", "peri", "i"]:
            print("WARNING: covariance labels are", cov.get("labels"), "- check the parameter order before use")
        (FRESH / f"{name}.json").write_text(json.dumps(obj, indent=1))
        print(f"wrote {FRESH / (name + '.json')}  (add the 'planet', 'esa', 'horizons' and 'analysis' blocks by hand; see data/objects/)")
    except (KeyError, TypeError, ValueError) as err:
        print(f"could not convert the SBDB record automatically ({err!r}); use the raw file")
    # --- JPL Horizons: osculating elements every 25 years, 1600-2500
    txt = get(HORIZONS, {"format": "text", "COMMAND": f"'DES={des};'", "OBJ_DATA": "NO", "EPHEM_TYPE": "ELEMENTS",
                         "CENTER": "'500@10'", "START_TIME": "1600-01-01", "STOP_TIME": "2500-01-01", "STEP_SIZE": "'25 y'",
                         "CSV_FORMAT": "YES", "OUT_UNITS": "AU-D"}).text
    (FRESH / f"{name}_horizons_raw.txt").write_text(txt)
    try:
        body = txt.split("$$SOE")[1].split("$$EOE")[0]
        out = FRESH / f"{name}_horizons.csv"
        with open(out, "w") as fh:
            fh.write(f"# JPL Horizons, {des}, heliocentric ecliptic J2000 osculating elements, AU-D\n# JDTDB, EC, IN, OM, W, MA, A\n")
            for line in body.strip().splitlines():
                c = [x.strip() for x in line.split(",")]
                # columns: JDTDB, date, EC, QR, IN, OM, W, Tp, N, MA, TA, A, AD, PR
                fh.write(", ".join([c[0], c[2], c[4], c[5], c[6], c[9], c[11]]) + "\n")
        print(f"wrote {out}")
    except (IndexError, ValueError) as err:
        print(f"could not parse the Horizons table ({err!r}); see the raw file")
    # --- ESA NEO Coordination Centre: independently computed orbit
    try:
        txt = get(NEOCC, {"file": des.replace(" ", "") + ".ke1"}).text
        (FRESH / f"{name}_esa_raw.txt").write_text(txt)
        print(f"saved ESA orbit file -> {FRESH / (name + '_esa_raw.txt')}  (lines KEP, MJD, MAG, RMS)")
    except Exception as err:  # noqa: BLE001
        print(f"ESA request failed ({err!r})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="what", required=True)
    z = sub.add_parser("zone")
    z.add_argument("zone", choices=sorted(ZONES))
    o = sub.add_parser("object")
    o.add_argument("designation", help='for example "2026 LO4"')
    args = ap.parse_args()
    if args.what == "zone":
        fetch_zone(args.zone)
    else:
        fetch_object(args.designation)


if __name__ == "__main__":
    main()
