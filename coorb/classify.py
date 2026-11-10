"""Reading co-orbital behaviour off the resonant angle phi(t) = lambda - lambda_planet.

A body in 1:1 resonance has phi that librates (swings back and forth) instead of
circulating. Where the swing sits tells the type:

    quasi-satellite (QS)   phi librates around 0        (stays beside the planet)
    Trojan / tadpole       phi librates around +60 (L4) or -60 (L5)
    horseshoe (HS)         phi librates around 180 with a gap at the planet
    compound               one swing encloses both 180 and 0 (e.g. horseshoe + QS)

Everything here works on a time series sampled at a fixed cadence, with t = 0 at the
element epoch. Angles are in degrees.
"""
import numpy as np

QS_CENTRE_LIMIT = 35.0     # deg; one-sided swings centred closer than this to the planet count as QS


def wrap180(x):
    return (np.asarray(x) + 180.0) % 360.0 - 180.0


def unwrap_from_epoch(t, phi):
    """Continuous (unwrapped) phi, shifted so it equals the wrapped value at t = 0."""
    i0 = int(np.argmin(np.abs(t)))
    pu = np.rad2deg(np.unwrap(np.deg2rad(np.asarray(phi, dtype=float))))
    return pu - pu[i0] + phi[i0]


def smooth(x, n):
    """Centred running mean over n samples (n odd), edges padded."""
    n = max(1, int(n)) | 1
    if n == 1:
        return np.asarray(x, dtype=float)
    return np.convolve(np.pad(x, n // 2, mode="edge"), np.ones(n) / n, mode="valid")


def samples_for(years, dt):
    """Odd number of samples spanning about `years` at cadence `dt`."""
    return max(1, int(round(years / dt))) | 1


def zigzag(t, y, prom):
    """Indices of the turning points of y(t), ignoring wiggles smaller than `prom`.

    A point counts as a turning point only if y moves by more than `prom` on both sides of
    it, so an extremum closer than that to the start of the series is not reported.
    """
    piv = []
    direction = 0
    i_min = i_max = i_ext = 0
    for i in range(1, len(y)):
        if direction == 0:                       # still looking for the first clear move
            if y[i] > y[i_max]:
                i_max = i
            if y[i] < y[i_min]:
                i_min = i
            if y[i] - y[i_min] > prom:
                direction = 1
                piv.append(i_min)
                i_ext = i
            elif y[i_max] - y[i] > prom:
                direction = -1
                piv.append(i_max)
                i_ext = i
        elif direction == 1:
            if y[i] > y[i_ext]:
                i_ext = i
            elif y[i_ext] - y[i] > prom:
                piv.append(i_ext)
                direction = -1
                i_ext = i
        else:
            if y[i] < y[i_ext]:
                i_ext = i
            elif y[i] - y[i_ext] > prom:
                piv.append(i_ext)
                direction = 1
                i_ext = i
    if piv and abs(y[piv[0]] - y[0]) <= prom:    # first pivot has no clear move before it
        piv = piv[1:]
    return piv


def turning_points(t, phi, smooth_n, prom=30.0):
    """Times and (unwrapped, smoothed) angles at which the drift of phi reverses."""
    ps = smooth(unwrap_from_epoch(t, phi), smooth_n)
    piv = zigzag(t, ps, prom)
    return t[piv], ps[piv], ps


def classify_now(t, phi, smooth_n=15, prom=8.0):
    """Label the state at t = 0 (used for screening many objects).

    The labels are deliberately coarse: they sort a few hundred objects so that the
    interesting ones can be looked at properly, and borderline cases can land either way.

    Returns a dict with `state` in {QS, L4, L5, HS, compound, passing, edge, unstable} and,
    for librating states, the two turning angles of the current swing (`lo`, `hi`) and the
    interval [`t_start`, `t_end`] over which every swing stays below a full turn.
    """
    t = np.asarray(t, dtype=float)
    phi = np.asarray(phi, dtype=float)
    ok = np.isfinite(phi)
    if ok.sum() < 0.9 * len(phi):          # ejected at some point: keep the stretch around now
        i0 = int(np.argmin(np.abs(t)))
        lo = i0
        while lo > 0 and ok[lo - 1]:
            lo -= 1
        hi = i0
        while hi < len(phi) - 1 and ok[hi + 1]:
            hi += 1
        t = t[lo:hi + 1]
        phi = phi[lo:hi + 1]
    else:
        phi = np.where(ok, phi, 0.0)
    if len(t) < 50:
        return dict(state="unstable")
    i0 = int(np.argmin(np.abs(t)))
    tp, pp, ps = turning_points(t, phi, smooth_n, prom)
    res = dict(phi0=float(phi[i0]), n_turn=len(tp))
    if len(tp) == 0:
        res.update(state="passing", span=0.0)
        return res
    before = [j for j in range(len(tp)) if tp[j] <= 0]
    after = [j for j in range(len(tp)) if tp[j] > 0]
    if not before or not after:            # a reversal on one side of now only
        j = after[0] if after else before[-1]
        edge = ps[0] if after else ps[-1]
        if abs(pp[j] - edge) > 360.0:
            res.update(state="passing", span=0.0)
        else:
            res.update(state="edge", span=float(t[-1] - t[0]))
        return res
    jb, ja = before[-1], after[0]
    leg = abs(pp[ja] - pp[jb])
    if leg >= 360.0:
        res.update(state="passing", span=0.0, leg=float(leg))
        return res
    j1 = jb
    while j1 > 0 and abs(pp[j1] - pp[j1 - 1]) < 360.0:
        j1 -= 1
    j2 = ja
    while j2 < len(tp) - 1 and abs(pp[j2 + 1] - pp[j2]) < 360.0:
        j2 += 1
    t_start = tp[j1] if j1 > 0 or abs(pp[j1] - ps[0]) >= 360.0 else t[0]
    t_end = tp[j2] if j2 < len(tp) - 1 or abs(pp[j2] - ps[-1]) >= 360.0 else t[-1]
    lo, hi = min(pp[jb], pp[ja]), max(pp[jb], pp[ja])

    def encloses(x):                       # is angle x (mod 360) inside the current swing?
        k = np.ceil((lo - x) / 360.0)
        return (x + 360.0 * k) <= hi

    enc0, enc180 = bool(encloses(0.0)), bool(encloses(180.0))
    centre = float(wrap180(0.5 * (lo + hi)))
    if enc0 and not enc180:
        state = "QS"
    elif enc180 and not enc0:
        state = "HS"
    elif not enc0 and not enc180:
        # a swing that stays on one side of the planet: Trojan-like if it is centred well away
        # from the planet, otherwise an off-centre quasi-satellite (35 deg threshold)
        state = "QS" if abs(centre) < QS_CENTRE_LIMIT else ("L4" if centre > 0 else "L5")
    else:
        state = "compound"
    res.update(state=state, leg=float(leg), lo=float(wrap180(lo)), hi=float(wrap180(hi)), centre=centre,
               t_start=float(t_start), t_end=float(t_end), span=float(t_end - t_start))
    return res


def bounded_interval(t, phi, smooth_n, prom=30.0):
    """Interval around t = 0 over which phi librates (every swing spans less than 360 deg).

    Returns dict(start, end, period, at_run_start, at_run_end). `start`/`end` equal the ends
    of the run when the libration is still going there (flagged by the booleans). NaNs if
    there is no reversal on both sides of t = 0.
    """
    tp, pp, ps = turning_points(t, phi, smooth_n, prom)
    nan = dict(start=np.nan, end=np.nan, period=np.nan, at_run_start=False, at_run_end=False)
    if len(tp) < 2:
        return nan
    legs = np.abs(np.diff(pp))
    jb = int(np.searchsorted(tp, 0)) - 1
    if jb < 0 or jb >= len(tp) - 1:
        return nan
    lo = jb
    while lo > 0 and legs[lo - 1] < 360:
        lo -= 1
    hi = jb
    while hi < len(legs) and legs[hi] < 360:
        hi += 1
    open_start = lo == 0 and abs(pp[0] - ps[0]) < 360
    open_end = hi == len(legs) and abs(ps[-1] - pp[-1]) < 360
    start = t[0] if open_start else tp[lo]
    end = t[-1] if open_end else tp[hi]
    period = tp[jb + 2] - tp[jb] if jb + 2 < len(tp) else np.nan
    return dict(start=float(start), end=float(end), period=float(period),
                at_run_start=bool(open_start), at_run_end=bool(open_end))


def pure_horseshoe_interval(t, phi, smooth_n):
    """Interval around t = 0 during which the (smoothed) phi never crosses the planet.

    phi is unwrapped so that the value at t = 0 lies in (0, 360); the interval is the longest
    stretch containing t = 0 in which it stays inside (0, 360). Also returns the closest
    angular approach to the planet on the leading and trailing side within that stretch.
    """
    t = np.asarray(t, dtype=float)
    i0 = int(np.argmin(np.abs(t)))
    pu = np.rad2deg(np.unwrap(np.deg2rad(np.asarray(phi, dtype=float))))
    pu = pu - pu[i0] + np.mod(phi[i0], 360.0)
    ps = smooth(pu, smooth_n)
    inside = (ps > 0) & (ps < 360)
    lo = i0
    while lo > 0 and inside[lo - 1]:
        lo -= 1
    hi = i0
    while hi < len(t) - 1 and inside[hi + 1]:
        hi += 1
    seg = ps[lo:hi + 1]
    return dict(start=float(t[lo]), end=float(t[hi]), at_run_start=bool(lo == 0), at_run_end=bool(hi == len(t) - 1),
                gap_ahead=float(seg.min()), gap_behind=float(360 - seg.max()))


def lane_switches(t, a, a_planet, smooth_n):
    """Times at which the smoothed (a - a_planet) changes sign: the asteroid changes lane."""
    s = np.sign(smooth(np.asarray(a) - np.asarray(a_planet), smooth_n))
    return np.asarray(t)[1:][np.diff(s) != 0]
