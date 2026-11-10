"""Planetary initial conditions and masses from the JPL DE421 ephemeris.

DE421 is read through the `jplephem` and `de421` packages from PyPI, so no download is
needed. It is used only to place the Sun and planets at the asteroid's epoch; from there
every body is integrated together. All vectors are barycentric, in the ecliptic-J2000 frame
that JPL uses for small-body orbital elements (IAU76 obliquity, 84381.448 arcsec).

Units: au, day. Masses are GM values in au^3/day^2, so the simulation uses G = 1.
"""
import numpy as np
import de421
from jplephem import Ephemeris

_EPH = Ephemeris(de421)
AU_KM = float(_EPH.AU)
_EPS = np.deg2rad(84381.448 / 3600.0)
_CE, _SE = np.cos(_EPS), np.sin(_EPS)

EMRAT = float(_EPH.EMRAT)     # Earth/Moon mass ratio
GM = {
    "sun": float(_EPH.GMS), "mercury": float(_EPH.GM1), "venus": float(_EPH.GM2),
    "emb": float(_EPH.GMB), "mars": float(_EPH.GM4), "jupiter": float(_EPH.GM5),
    "saturn": float(_EPH.GM6), "uranus": float(_EPH.GM7), "neptune": float(_EPH.GM8),
}
GM["earth"] = GM["emb"] * EMRAT / (1.0 + EMRAT)
GM["moon"] = GM["emb"] / (1.0 + EMRAT)

# DE421 time coverage (Julian dates)
JD_MIN, JD_MAX = float(_EPH.jalpha), float(_EPH.jomega)


def equatorial_to_ecliptic(v):
    x, y, z = v
    return np.array([x, _CE * y + _SE * z, -_SE * y + _CE * z])


def bary_state(name, jd):
    """Barycentric ecliptic-J2000 position [au] and velocity [au/day] of a body.

    `name` is one of: sun, mercury, venus, earth, moon, emb (Earth-Moon barycentre), mars,
    jupiter, saturn, uranus, neptune. Outer planets are system barycentres.
    """
    if not (JD_MIN <= jd <= JD_MAX):
        raise ValueError(f"JD {jd} is outside DE421 coverage ({JD_MIN}..{JD_MAX})")
    if name in ("earth", "moon"):
        pe, ve = _EPH.position_and_velocity("earthmoon", jd)
        pm, vm = _EPH.position_and_velocity("moon", jd)      # geocentric Moon
        pe, ve, pm, vm = (np.asarray(a).ravel() for a in (pe, ve, pm, vm))
        p_earth = pe - pm / (1.0 + EMRAT)
        v_earth = ve - vm / (1.0 + EMRAT)
        if name == "earth":
            p, v = p_earth, v_earth
        else:
            p, v = p_earth + pm, v_earth + vm
    else:
        key = "earthmoon" if name == "emb" else name
        p, v = _EPH.position_and_velocity(key, jd)
        p, v = np.asarray(p).ravel(), np.asarray(v).ravel()
    return equatorial_to_ecliptic(p / AU_KM), equatorial_to_ecliptic(v / AU_KM)
