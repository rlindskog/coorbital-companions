"""Fast checks of the building blocks. Run with:  python -m pytest -q"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from coorb import classify, ephem, nbody, orbits  # noqa: E402


def test_elements_round_trip():
    el = (1.3, 0.2, 10.0, 50.0, 120.0, 250.0)
    pos, vel = nbody.kepler_to_state(*el)
    out = nbody.full_elements(pos, vel)
    got = (out["a"], out["e"], out["i"], out["node"], out["peri"], out["M"])
    assert np.allclose(got, el, rtol=0, atol=1e-9)
    a, e, inc, lam = nbody.elements(ephem.GM["sun"], pos, vel)
    assert abs(np.degrees(lam) - (50.0 + 120.0 + 250.0) % 360) < 1e-9


def test_ephemeris_is_sane():
    jd = 2461200.5
    ps, vs = ephem.bary_state("sun", jd)
    pe, ve = ephem.bary_state("earth", jd)
    assert 0.98 < np.linalg.norm(pe - ps) < 1.02
    pv, vv = ephem.bary_state("venus", jd)
    a = nbody.full_elements(pv - ps, vv - vs)["a"]
    assert abs(a - 0.7233) < 0.001
    assert abs(pe[2] - ps[2]) < 1e-3          # Earth stays in the ecliptic plane of the frame


def test_model_reproduces_jpl_propagation():
    """Element epoch -> covariance epoch must land on JPL's own elements (see check_propagation.py)."""
    for name in ("2026_LO4", "2025_ON9"):
        obj = orbits.load_object(name)
        jpl, cov = obj["jpl"], obj["jpl"]["covariance"]
        model = nbody.Model(jpl["epoch_jd_tdb"], [orbits.element_tuple(jpl["elements"])], gr=True, split_moon=True)
        model.advance(cov["epoch_jd_tdb"] - jpl["epoch_jd_tdb"])
        el = model.helio_elements(0)
        mine = dict(e=el["e"], q=el["q"], node=el["node"], peri=el["peri"], i=el["i"])
        for k, lab in enumerate(cov["labels"]):
            if lab == "tp":
                continue
            sigma = cov["matrix"][k][k] ** 0.5
            assert abs(mine[lab] - cov["nominal"][k]) < 0.01 * sigma, (name, lab)


def test_clone_sampling_matches_covariance():
    obj = orbits.load_object("2026_LO4")
    epoch, asts = orbits.sample_clones(obj, 4000, seed=3)
    assert epoch == obj["jpl"]["covariance"]["epoch_jd_tdb"]
    a = np.array([x[0] for x in asts[1:]])
    assert abs(a.std() / obj["jpl"]["sigma"]["a_au"] - 1) < 0.08
    nominal = orbits.kepler_from_cov_params(obj["jpl"]["covariance"]["nominal"], epoch)
    assert asts[0] == nominal


def test_classifier_on_synthetic_signals():
    t = np.arange(-3000.0, 3000.0, 0.5)
    w = lambda x: classify.wrap180(x)
    n = 7
    horseshoe = w(180 + 160 * np.sin(2 * np.pi * t / 400 + 1.0))
    quasi = 15 * np.sin(2 * np.pi * t / 90)
    trojan = 60 + 25 * np.sin(2 * np.pi * t / 220)
    passing = w(0.9 * t)
    assert classify.classify_now(t, horseshoe, smooth_n=n)["state"] == "HS"
    assert classify.classify_now(t, quasi, smooth_n=n)["state"] == "QS"
    assert classify.classify_now(t, trojan, smooth_n=n)["state"] == "L4"
    assert classify.classify_now(t, w(-trojan), smooth_n=n)["state"] == "L5"
    assert classify.classify_now(t, passing, smooth_n=n)["state"] == "passing"
    off_centre = -14 + 13 * np.sin(2 * np.pi * t / 150)          # beside the planet, never crossing it
    assert classify.classify_now(t, off_centre, smooth_n=n)["state"] == "QS"
    hs = classify.pure_horseshoe_interval(t, horseshoe, n)
    assert hs["at_run_start"] and hs["at_run_end"] and 15 < hs["gap_ahead"] < 25
    lib = classify.bounded_interval(t, horseshoe, n)
    assert lib["at_run_start"] and lib["at_run_end"] and abs(lib["period"] - 400) < 5
    # a libration that starts 1000 years ago, after circulating
    mixed = np.where(t < -1000, w(1.2 * (t + 1000) + 180 + 160 * np.sin(1.0 - 2 * np.pi * 1000 / 400)), horseshoe)
    lib = classify.bounded_interval(t, mixed, n)
    assert not lib["at_run_start"] and -1300 < lib["start"] < -700


def test_horizons_agreement_near_epoch():
    """The model should sit on JPL Horizons' own trajectory (full comparison: compare_horizons.py)."""
    obj = orbits.load_object("2026_LO4")
    H = orbits.load_horizons("2026_LO4")
    row = H[np.argmin(np.abs(H[:, 0] - 2469807.5))]          # 2050-Jan-01
    ep0 = obj["jpl"]["epoch_jd_tdb"]
    model = nbody.Model(ep0, [orbits.element_tuple(obj["jpl"]["elements"])], gr=True, split_moon=True)
    model.advance(row[0] - ep0)
    el = model.helio_elements(0)
    lam_model = (el["node"] + el["peri"] + el["M"]) % 360
    lam_h = (row[3] + row[4] + row[5]) % 360
    assert abs((lam_model - lam_h + 180) % 360 - 180) < 1e-4
    assert abs(el["a"] - row[6]) < 1e-8


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
