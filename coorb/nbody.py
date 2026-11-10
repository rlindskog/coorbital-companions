"""N-body model: Sun and planets as massive bodies, asteroids as massless test particles.

Integrator: REBOUND's IAS15 (adaptive, 15th order). Optional extras through REBOUNDx:
general relativity (post-Newtonian, Sun only) and a slow drift in semi-major axis that
mimics the Yarkovsky effect.

Conventions
-----------
* Units: au, day; masses are GM values (G = 1).
* Asteroid elements are heliocentric, ecliptic J2000, osculating at a TDB epoch, as
  published by the JPL Small-Body Database: (a [au], e, i, node, peri, M [deg]).
* The resonant angle reported everywhere is the difference of mean longitudes,
  phi = lambda_asteroid - lambda_planet, in degrees in (-180, 180].
"""
import numpy as np
import rebound

from .ephem import GM, bary_state

INNER = ["mercury", "venus", "emb", "mars"]
OUTER = ["jupiter", "saturn", "uranus", "neptune"]
C_AU_PER_DAY = 173.14463267424034      # speed of light
DAYS_PER_YEAR = 365.25


# --------------------------------------------------------------------------- set-up
def build_simulation(jd, bodies=None, split_moon=False):
    """Sun + planets at Julian date `jd` (TDB). Returns (sim, names).

    `bodies` defaults to all eight planets, with Earth and Moon merged at their barycentre
    ("emb"). `split_moon=True` carries Earth and Moon as separate bodies instead. Planets
    left out of `bodies` have their mass added to the Sun.
    """
    if bodies is None:
        bodies = INNER + OUTER
    bodies = list(bodies)
    if split_moon and "emb" in bodies:
        k = bodies.index("emb")
        bodies[k:k + 1] = ["earth", "moon"]
    sim = rebound.Simulation()
    sim.G = 1.0
    msun = GM["sun"]
    for b in INNER + OUTER:
        if b not in bodies and not (b == "emb" and split_moon):
            msun += GM[b]
    p, v = bary_state("sun", jd)
    sim.add(m=msun, x=p[0], y=p[1], z=p[2], vx=v[0], vy=v[1], vz=v[2])
    for b in bodies:
        p, v = bary_state(b, jd)
        sim.add(m=GM[b], x=p[0], y=p[1], z=p[2], vx=v[0], vy=v[1], vz=v[2])
    sim.N_active = sim.N
    sim.testparticle_type = 0
    sim.integrator = "ias15"
    return sim, ["sun"] + bodies


def kepler_to_state(a, e, inc, node, peri, M, mu=None):
    """Heliocentric position and velocity from osculating elements (angles in degrees).

    Uses mu = GM_sun, the convention behind JPL's heliocentric elements for a massless body.
    """
    mu = GM["sun"] if mu is None else mu
    inc, node, peri, M = np.deg2rad([inc, node, peri, M])
    E = M
    for _ in range(100):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= dE
        if abs(dE) < 1e-15:
            break
    xp = a * (np.cos(E) - e)
    yp = a * np.sqrt(1 - e * e) * np.sin(E)
    r = a * (1 - e * np.cos(E))
    n = np.sqrt(mu / a ** 3)
    vxp = -a * n * np.sin(E) * a / r
    vyp = a * n * np.sqrt(1 - e * e) * np.cos(E) * a / r
    cO, sO, cw, sw, ci, si = np.cos(node), np.sin(node), np.cos(peri), np.sin(peri), np.cos(inc), np.sin(inc)
    P = np.array([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si])
    Q = np.array([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si])
    return xp * P + yp * Q, vxp * P + vyp * Q


def add_asteroid(sim, a, e, inc, node, peri, M):
    """Add a massless particle from heliocentric elements (the Sun is particle 0)."""
    s = sim.particles[0]
    pos, vel = kepler_to_state(a, e, inc, node, peri, M)
    sim.add(m=0.0, x=s.x + pos[0], y=s.y + pos[1], z=s.z + pos[2],
            vx=s.vx + vel[0], vy=s.vy + vel[1], vz=s.vz + vel[2])


# --------------------------------------------------------------------------- elements
def elements(mu, r, v):
    """Vectorised a, e, inc [rad], mean longitude [rad] from relative states (..., 3)."""
    r = np.asarray(r)
    v = np.asarray(v)
    rn = np.linalg.norm(r, axis=-1)
    v2 = np.sum(v * v, axis=-1)
    h = np.cross(r, v)
    hn = np.linalg.norm(h, axis=-1)
    with np.errstate(all="ignore"):
        a = 1.0 / (2.0 / rn - v2 / mu)
        rv = np.sum(r * v, axis=-1)
        ecosf = hn * hn / (mu * rn) - 1.0
        esinf = rv * hn / (mu * rn)
        e = np.hypot(ecosf, esinf)
        f = np.arctan2(esinf, ecosf)
        inc = np.arccos(np.clip(h[..., 2] / hn, -1, 1))
        nx, ny = -h[..., 1], h[..., 0]
        node = np.arctan2(ny, nx)
        cu = r[..., 0] * nx + r[..., 1] * ny
        hxn = np.stack([-h[..., 2] * ny, h[..., 2] * nx, h[..., 0] * ny - h[..., 1] * nx], axis=-1)
        su = np.sum(r * hxn, axis=-1) / hn
        u = np.arctan2(su, cu)                       # argument of latitude
        E = 2.0 * np.arctan2(np.sqrt(np.clip(1 - e, 0, None)) * np.sin(f / 2),
                             np.sqrt(1 + e) * np.cos(f / 2))
        M = E - e * np.sin(E)
        lam = node + u - f + M
        bad = ~(e < 1.0)
        a = np.where(bad, np.nan, a)
        lam = np.where(bad, np.nan, lam)
    return a, e, inc, np.mod(lam, 2 * np.pi)


def full_elements(r, v, mu=None):
    """All six osculating elements of one body (angles in degrees) plus q and n [deg/day]."""
    mu = GM["sun"] if mu is None else mu
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    h = np.cross(r, v)
    hn = np.linalg.norm(h)
    rn = np.linalg.norm(r)
    a = 1 / (2 / rn - v @ v / mu)
    ev = np.cross(v, h) / mu - r / rn
    e = np.linalg.norm(ev)
    inc = np.degrees(np.arccos(h[2] / hn))
    n = np.array([-h[1], h[0], 0.0])
    node = np.degrees(np.arctan2(n[1], n[0])) % 360
    peri = np.degrees(np.arctan2(np.dot(np.cross(n, ev), h) / hn, np.dot(n, ev))) % 360
    f = np.arctan2(np.dot(np.cross(ev, r), h) / hn, np.dot(ev, r))
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(f / 2), np.sqrt(1 + e) * np.cos(f / 2))
    M = np.degrees(E - e * np.sin(E)) % 360
    return dict(a=a, e=e, i=inc, node=node, peri=peri, M=M, q=a * (1 - e),
                n=np.degrees(np.sqrt(mu / a ** 3)))


# --------------------------------------------------------------------------- model
class Model:
    """Sun + planets + test particles, ready to step to arbitrary times.

    Parameters
    ----------
    jd : float
        Epoch (JD, TDB) of the asteroid elements; time zero of the simulation.
    asteroids : sequence of (a, e, i, node, peri, M)
    gr : bool
        Add the Sun's post-Newtonian correction (REBOUNDx "gr").
    split_moon : bool
        Carry Earth and Moon separately rather than as one body at their barycentre.
    drift_au_per_myr : sequence of float or None
        One da/dt value per asteroid (au per million years); 0 or None means no drift.
    """

    def __init__(self, jd, asteroids, *, gr=True, split_moon=False, bodies=None, drift_au_per_myr=None):
        self.jd0 = float(jd)
        self.sim, self.names = build_simulation(jd, bodies=bodies, split_moon=split_moon)
        self.n_massive = self.sim.N
        asteroids = [tuple(el) for el in asteroids]
        for el in asteroids:
            add_asteroid(self.sim, *el)
        self.n_asteroids = len(asteroids)
        self.sim.N_active = self.n_massive
        self._rebx = None                      # keep a reference alive for the run
        has_drift = drift_au_per_myr is not None and any(r for r in drift_au_per_myr)
        if gr or has_drift:
            import reboundx
            self._rebx = reboundx.Extras(self.sim)
            if gr:
                g = self._rebx.load_force("gr")
                g.params["c"] = C_AU_PER_DAY
                self._rebx.add_force(g)
            if has_drift:
                mof = self._rebx.load_force("modify_orbits_forces")
                self._rebx.add_force(mof)
                mof.params["coordinates"] = reboundx.coordinates["PARTICLE"]
                self.sim.particles[0].params["primary"] = 1
                for k, rate in enumerate(drift_au_per_myr):
                    if rate:
                        a0 = asteroids[k][0]
                        self.sim.particles[self.n_massive + k].params["tau_a"] = a0 / (rate / 1e6 / DAYS_PER_YEAR)
        self.sim.move_to_com()
        self._buf = np.zeros((self.sim.N, 6))

    def advance(self, t_days):
        """Integrate to `t_days` from the epoch (negative = past); returns the state array."""
        self.sim.integrate(t_days, exact_finish_time=1)
        self.sim.serialize_particle_data(xyzvxvyvz=self._buf)
        return self._buf

    def index(self, name):
        return self.names.index(name)

    def planet_state(self, planet):
        """Barycentric state of a planet; Earth is the Earth-Moon barycentre when split."""
        buf = self._buf
        if planet in ("earth", "emb") and "earth" in self.names:
            me, mm = GM["earth"], GM["moon"]
            return (me * buf[self.index("earth")] + mm * buf[self.index("moon")]) / (me + mm)
        return buf[self.index(planet)]

    def asteroid_states(self):
        return self._buf[self.n_massive:]

    def helio_elements(self, k=0):
        """Six osculating heliocentric elements of asteroid k at the current time."""
        sun = self._buf[0]
        s = self._buf[self.n_massive + k]
        return full_elements(s[:3] - sun[:3], s[3:] - sun[3:])


# --------------------------------------------------------------------------- runs
def integrate(jd, asteroids, planet, years, dt_out_years=0.2, direction=+1, **model_kwargs):
    """Integrate for `years` in one direction, sampling every `dt_out_years`.

    Returns a dict of arrays: t [years from epoch], a [au], e, inc [rad], phi [deg], all of
    shape (n_times, n_asteroids); a_pl [au] (planet's semi-major axis) and d [au] (distance
    from each asteroid to the planet).
    """
    model = Model(jd, asteroids, **model_kwargs)
    nast = model.n_asteroids
    nt = int(round(years / dt_out_years)) + 1
    t_out = direction * np.arange(nt) * dt_out_years * DAYS_PER_YEAR
    A = np.full((nt, nast), np.nan)
    EC = np.full((nt, nast), np.nan)
    IN = np.full((nt, nast), np.nan)
    PHI = np.full((nt, nast), np.nan)
    D = np.full((nt, nast), np.nan)
    APL = np.full(nt, np.nan)
    mu = GM["sun"]
    for k, t in enumerate(t_out):
        buf = model.advance(t)
        sun = buf[0]
        pl = model.planet_state(planet)
        ap, ep, incp, lamp = elements(mu, pl[:3] - sun[:3], pl[3:] - sun[3:])
        rel = buf[model.n_massive:] - sun
        a, e, inc, lam = elements(mu, rel[:, :3], rel[:, 3:])
        A[k] = a
        EC[k] = e
        IN[k] = inc
        PHI[k] = np.rad2deg(np.mod(lam - lamp + np.pi, 2 * np.pi) - np.pi)
        APL[k] = ap
        D[k] = np.linalg.norm(buf[model.n_massive:, :3] - pl[:3], axis=1)
    return dict(t=t_out / DAYS_PER_YEAR, a=A, e=EC, inc=IN, phi=PHI, a_pl=APL, d=D)


def join_directions(backward, forward):
    """Concatenate a backward and a forward run into one series in increasing time."""
    return {k: np.concatenate([backward[k][::-1][:-1], forward[k]], axis=0) for k in forward}


def integrate_both(jd, asteroids, planet, years, **kwargs):
    """Backward and forward from the epoch, in sequence, joined in increasing time."""
    fw = integrate(jd, asteroids, planet, years, direction=+1, **kwargs)
    bw = integrate(jd, asteroids, planet, years, direction=-1, **kwargs)
    return join_directions(bw, fw)


def _one_direction(args):
    jd, asteroids, planet, years, direction, kwargs = args
    return integrate(jd, asteroids, planet, years, direction=direction, **kwargs)


def integrate_both_parallel(jd, asteroids, planet, years, **kwargs):
    """As `integrate_both`, with the two directions on two processes."""
    from multiprocessing import Pool
    with Pool(2) as pool:
        fw, bw = pool.map(_one_direction, [(jd, asteroids, planet, years, +1, kwargs),
                                           (jd, asteroids, planet, years, -1, kwargs)])
    return join_directions(bw, fw)
