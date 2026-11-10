"""Draw the summary figure for each object from the files in results/.

Needs results/<object>/nominal_series.csv (scripts/nominal.py) and
results/<object>/horizons_comparison.csv (scripts/compare_horizons.py).

    python scripts/make_figures.py                # both objects
    python scripts/make_figures.py 2025_ON9
"""
import _bootstrap
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch

from coorb import classify, orbits
from coorb.paths import FIGURES, RESULTS, year_of_jd

SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8983", "#e4e3de"
BLUE, ORANGE = "#2a78d6", "#eb6834"

FIG = {
    "2026_LO4": dict(
        title="Asteroid 2026 LO4 is locked in step with Venus",
        subtitle="Every 330 years or so it drifts all the way around Venus’s orbit and back again.",
        planet="Venus", exaggeration=45.0, a_planet=0.7233, own_period="224-day",
        arrows=[(0, 0.23), (0, 0.44), (0, 0.65), (1, 0.25), (1, 0.50), (1, 0.74)], arrow_step=8,
        planet_label=(1.36, -0.19), now_label=(0.56, -0.36, "right"),
        turn_a=dict(text="other turning point:\n{y1} and {y3}", at=1, xy=(1.30, 0.80)),
        turn_b=dict(text="swings past Venus,\nturns back in {y2}", at=2, xy=(1.36, 0.20)),
        inner_label=(-0.40, 0.42), x_range=(1600, 2500), x_ticks=[1600, 1800, 2000, 2200, 2400], x_tick_labels=None,
        y_lim=(-12, 408), y_tick_labels=["at Venus", "90° ahead", "far side\nof the Sun", "270° ahead", "at Venus\n(lapped)"],
        band=(2014, 2071, "passing Venus\nnow, 2014–2071"), horizons_every=1, panel_b_title="WHERE IT SITS ALONG VENUS’S ORBIT, 1600–2500"),
    "2025_ON9": dict(
        title="2025 ON9 rides a horseshoe around Mars",
        subtitle="It swings from just ahead of Mars to just behind it and back every 1,750 years or so.",
        planet="Mars", exaggeration=60.0, a_planet=1.5237, own_period="1.9-year",
        arrows=[(0, 0.28), (0, 0.52), (0, 0.76), (1, 0.28), (1, 0.52), (1, 0.76)], arrow_step=30,
        planet_label=(1.42, 0.0), now_label=(0.05, 0.50, "center"),
        turn_a=dict(text="turned back here\nin {y1} and {y3}", at=1, xy=(1.42, 0.62)),
        turn_b=dict(text="turns back here\nin {y2}", at=2, xy=(1.42, -0.62)),
        inner_label=(-0.36, -0.42), x_range=(0, 4500), x_ticks=[0, 1000, 2000, 3000, 4000],
        x_tick_labels=["AD 1", "1000", "2000", "3000", "4000"],
        y_lim=(-12, 372), y_tick_labels=["at Mars", "90° ahead", "far side\nof the Sun", "90° behind", "at Mars"],
        band=None, horizons_every=2, panel_b_title="WHERE IT SITS ALONG MARS’S ORBIT, AD 1–4500"),
}


def draw(name):
    cfg = FIG[name]
    obj = orbits.load_object(name)
    an = obj["analysis"]
    S = np.genfromtxt(RESULTS / name / "nominal_series.csv", delimiter=",", names=True)
    H = np.genfromtxt(RESULTS / name / "horizons_comparison.csv", delimiter=",", names=True)
    yr, phi, a, a_pl = S["year"], S["phi_deg"], S["a_au"], S["a_planet_au"]
    now = year_of_jd(float(obj["jpl"]["epoch_jd_tdb"]))
    i0 = int(np.argmin(np.abs(yr - now)))
    pu = np.rad2deg(np.unwrap(np.deg2rad(phi)))
    pu = pu - pu[i0] + np.mod(phi[i0], 360.0)             # degrees ahead of the planet, continuous
    n = classify.samples_for(an["figure_smooth_years"], float(yr[1] - yr[0]))
    pus, da = classify.smooth(pu, n), classify.smooth(a - a_pl, n)
    piv = classify.zigzag(yr, pus, 30.0)
    jprev = max(j for j in range(len(piv)) if yr[piv[j]] <= now)
    j1, j2, j3 = piv[jprev], piv[jprev + 1], piv[jprev + 2]
    Y1, Y2, Y3 = (int(round(yr[j])) for j in (j1, j2, j3))
    K, AP, planet = cfg["exaggeration"], cfg["a_planet"], cfg["planet"]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.edgecolor": GRID, "axes.linewidth": 0.8,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.labelcolor": INK2})
    fig = plt.figure(figsize=(6.5, 9.6), dpi=200, facecolor=SURF)
    fig.text(0.06, 0.972, cfg["title"], fontsize=15.5, fontweight="bold", va="top")
    fig.text(0.06, 0.944, cfg["subtitle"], fontsize=10.2, color=INK2, va="top")

    # ---- panel A: one full cycle in a frame that turns with the planet
    ax = fig.add_axes([0.02, 0.475, 0.96, 0.425])
    ax.set_facecolor(SURF)
    ax.set_aspect("equal")
    ax.axis("off")
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color=GRID, lw=1.0, zorder=1)

    def xy(j):
        r = 1 + K * da[j] / AP
        return r * np.cos(np.deg2rad(pus[j])), r * np.sin(np.deg2rad(pus[j]))

    seg = np.arange(j1, j3 + 1)
    r = 1 + K * da[seg] / AP
    ax.plot(r * np.cos(np.deg2rad(pus[seg])), r * np.sin(np.deg2rad(pus[seg])), color=BLUE, lw=2.0,
            solid_capstyle="round", solid_joinstyle="round", zorder=3)
    for half, frac in cfg["arrows"]:
        ya, yb = (Y1, Y2) if half == 0 else (Y2, Y3)
        j = int(np.argmin(np.abs(yr - (ya + frac * (yb - ya)))))
        ax.add_patch(FancyArrowPatch(xy(j), xy(j + cfg["arrow_step"]), arrowstyle="-|>", mutation_scale=17, color=BLUE, lw=0, zorder=4))
    ax.add_patch(Circle((0, 0), 0.045, color=INK2, zorder=5))
    ax.text(0, -0.11, "Sun", ha="center", va="top", fontsize=10, color=INK2)
    ax.add_patch(Circle((1, 0), 0.036, facecolor=INK, edgecolor=SURF, lw=2, zorder=6))
    line = dict(arrowstyle="-", color=MUTED, lw=0.8)
    ax.annotate(planet, xy=(1.0, 0.0), xytext=cfg["planet_label"], fontsize=10.5, fontweight="bold", ha="left", va="center",
                arrowprops=dict(shrinkA=2, shrinkB=5, **line))
    ax.add_patch(Circle(xy(i0), 0.036, facecolor=BLUE, edgecolor=SURF, lw=2, zorder=7))
    nx, ny, nha = cfg["now_label"]
    ax.annotate(f"{obj['designation']} today", xy=xy(i0), xytext=(nx, ny), fontsize=10.5, fontweight="bold", ha=nha, va="center",
                arrowprops=dict(shrinkA=2, shrinkB=6, **line))
    for key in ("turn_a", "turn_b"):
        t = cfg[key]
        j = j1 if t["at"] == 1 else j2
        ax.annotate(t["text"].format(y1=Y1, y2=Y2, y3=Y3), xy=(np.cos(np.deg2rad(pus[j])), np.sin(np.deg2rad(pus[j]))), xytext=t["xy"],
                    fontsize=10, color=INK2, ha="left", va="center", arrowprops=dict(shrinkA=2, shrinkB=3, **line))
    ax.text(*cfg["inner_label"], f"inner lane:\ngains on {planet}", fontsize=9.5, color=INK2, ha="center", va="center")
    ax.text(-1.56, 0.86, "outer lane:\nfalls behind", fontsize=9.5, color=INK2, ha="center", va="center")
    ax.annotate("", xy=(-0.90, 0.68), xytext=(-1.38, 0.80), arrowprops=line)
    ax.set_xlim(-1.95, 2.35)
    ax.set_ylim(-1.38, 1.38)
    fig.text(0.06, 0.905, f"ONE FULL CYCLE ({Y1}–{Y3}), AS SEEN RIDING ALONG WITH {planet.upper()}", fontsize=8.5, color=MUTED, va="top")
    fig.text(0.06, 0.478, f"Lane spacing exaggerated {K:.0f}× so the two lanes can be seen. Each point is the asteroid’s average\n"
             f"position; its own {cfg['own_period']} loop around the Sun is smoothed out.", fontsize=8.2, color=MUTED, va="top", linespacing=1.45)

    # ---- panel B: degrees ahead of the planet against time, with JPL Horizons on top
    bx = fig.add_axes([0.20, 0.085, 0.74, 0.275])
    bx.set_facecolor(SURF)
    x0, x1 = cfg["x_range"]
    m = (yr >= x0) & (yr <= x1)
    if cfg["band"]:
        b0, b1, btxt = cfg["band"]
        bx.axvspan(b0, b1, color=BLUE, alpha=0.10, lw=0, zorder=0)
        bx.text(b1 + 11, cfg["y_lim"][1] - 8, btxt, fontsize=9.2, color=INK2, ha="left", va="top", linespacing=1.3)
    else:
        bx.axvline(now, color=MUTED, lw=0.8, zorder=2)
        bx.text(now + 40, 12, "today", fontsize=9.2, color=INK2, ha="left", va="bottom")
    for yy in (0, 90, 180, 270, 360):
        bx.axhline(yy, color=GRID, lw=0.8, zorder=1)
    curve = pu if cfg["band"] else pus
    bx.plot(yr[m], curve[m], color=BLUE, lw=2.0, zorder=3, solid_joinstyle="round", label="This simulation")
    ref = np.interp(H["year"], yr, pu)
    hp = H["deg_ahead_of_planet"] + 360.0 * np.round((ref - H["deg_ahead_of_planet"]) / 360.0)
    e = cfg["horizons_every"]
    bx.plot(H["year"][::e], hp[::e], "o", ms=6.2, mfc=ORANGE, mec=SURF, mew=1.4, zorder=4, label="NASA JPL Horizons’ own calculation")
    bx.set_xlim(x0, x1)
    bx.set_ylim(*cfg["y_lim"])
    bx.set_yticks([0, 90, 180, 270, 360])
    bx.set_yticklabels(cfg["y_tick_labels"], fontsize=9.2)
    bx.set_xticks(cfg["x_ticks"])
    if cfg["x_tick_labels"]:
        bx.set_xticklabels(cfg["x_tick_labels"])
    bx.tick_params(length=0, labelsize=9.5)
    for s in ("top", "right", "left"):
        bx.spines[s].set_visible(False)
    bx.legend(loc="lower left", bbox_to_anchor=(-0.215, 1.035), frameon=False, fontsize=9.3, handlelength=1.5, borderaxespad=0,
              ncol=2, columnspacing=1.6, handletextpad=0.5, labelcolor=INK2)
    fig.text(0.06, 0.418, cfg["panel_b_title"], fontsize=8.5, color=MUTED, va="top")
    j = obj["jpl"]
    first, last = j["first_obs"][:7], j["last_obs"][:7]
    fig.text(0.06, 0.014, f"Orbit: NASA/JPL Small-Body Database, solution #{j['solution_id']} ({j['n_obs']} observations, {first} to {last}).\n"
             "Model: Sun, eight planets and general relativity. ESA’s independent orbit gives the same result.",
             fontsize=7.8, color=MUTED, va="bottom", linespacing=1.45)
    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / f"{name}.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    plt.close(fig)
    print(f"wrote {_bootstrap.show(out)}  (cycle {Y1} -> {Y2} -> {Y3})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("objects", nargs="*", default=sorted(FIG))
    args = ap.parse_args()
    for name in args.objects:
        draw(name)


if __name__ == "__main__":
    main()
