# Two asteroids in step with a planet: 2026 LO4 (Venus) and 2025 ON9 (Mars)

*Analysis of 10 November 2026. Unreviewed; see [Limitations](#8-limitations).*

## 1. Summary

A screen of the JPL Small-Body Database for asteroids that share a planet's orbital period
turned up two objects whose resonant behaviour does not appear in the published catalogues.

- **2026 LO4** is in 1:1 mean-motion resonance with **Venus**. It follows a compound
  horseshoe / quasi-satellite path with a cycle of about 330 years, and since 2014 it has
  been in the quasi-satellite part of that cycle, looping beside Venus until about 2071. Its
  eccentricity, 0.19, is about half that of any Venus co-orbital in the published lists, and
  it is the only one whose orbit stays entirely inside Earth's. It is large for the class,
  roughly 0.4–0.9 km across.
- **2025 ON9** is in 1:1 resonance with **Mars**, on a regular horseshoe orbit with a cycle
  of about 1,750 years. It is roughly 55–125 m across.

For both, the simulation used here reproduces NASA JPL Horizons' own integration at 37 dates
between 1600 and 2500, ESA's independently computed orbit gives the same behaviour, and the
conclusion holds for every orbit drawn from the JPL covariance matrix.

What is new, if anything, is the identification of these dynamical states. The asteroids
themselves were found by sky surveys, and whether someone has already noticed how they move
could not be ruled out (section 7).

## 2. Background

An asteroid whose orbital period equals a planet's is *co-orbital* with it. Seen from a frame
that rotates with the planet, such an object does not drift steadily past: the planet's
gravity makes it swing back and forth. Three basic patterns exist, plus mixtures.

| Type | What it does | Resonant angle φ |
|---|---|---|
| Trojan (tadpole) | stays near a point 60° ahead of or behind the planet | swings around +60° or −60° |
| Horseshoe | drifts round almost the whole orbit, turning back each time it nears the planet | swings around 180°, never reaching 0° |
| Quasi-satellite | stays beside the planet and appears to circle it, without being bound to it | swings around 0° |
| Compound | one swing combines two of the above, for example horseshoe and quasi-satellite | one swing encloses both 180° and 0° |

Here φ is the asteroid's mean longitude minus the planet's. Positive φ means the asteroid is
ahead of the planet along the orbit.

A horseshoe works through a lane change. On the inner lane (slightly smaller orbit) the
asteroid moves faster than the planet and gains on it. As it closes in from behind, the
planet's pull lifts it to the outer lane, where it is slower, and it falls away again. The
reverse happens at the other end.

## 3. Data and method

### 3.1 Candidates

The JPL Small-Body Database was queried on 10 November 2026 for every object with a
semi-major axis inside a planet's co-orbital zone, taken as the planet's semi-major axis
± (m_planet / 3 M_sun)^(1/3) of itself, and an orbit condition code of 5 or better.

| Planet | Semi-major axis range (au) | Objects |
|---|---|---|
| Venus | 0.7165 – 0.7301 | 20 |
| Earth | 0.990 – 1.010 | 136 |
| Mars | 1.5164 – 1.5310 | 134 |

The lists are in `data/zones/`.

### 3.2 Dynamical model

- Sun and eight planets as point masses; asteroids as massless particles.
- Planet positions, velocities and masses at the asteroid's epoch from JPL DE421. From there
  everything is integrated together.
- REBOUND's IAS15 integrator (adaptive step, 15th order).
- Screening runs are Newtonian with Earth and Moon merged. All runs on the two objects add
  the Sun's relativistic correction (REBOUNDx `gr`); the comparisons with JPL also carry
  Earth and Moon separately.
- Not modelled: perturbations by other asteroids, the Sun's oblateness, and, except in the
  test of section 5.5, thermal (Yarkovsky) forces.

### 3.3 Reading the simulations

For each object the series φ(t) is unwrapped, lightly smoothed, and its turning points are
located. Three quantities are used throughout.

- **Librating interval:** the stretch of time around the present in which every swing of φ
  between consecutive turning points is smaller than 360°. Outside it the asteroid circulates
  past the planet and is not in resonance.
- **Pure-horseshoe interval:** the stretch around the present in which φ never crosses 0°.
- **Lane switches:** the times when the smoothed semi-major axis crosses the planet's.

The screening labels (`coorb/classify.py`) are deliberately coarse. They sort a few hundred
objects so the interesting ones can be examined properly.

### 3.4 Checks applied to each object

1. **Internal consistency of JPL's orbit.** JPL gives the elements at a standard epoch and
   again at the covariance epoch. Integrating one to the other tests the model and the data.
2. **JPL Horizons.** Horizons' osculating elements every 25 years from 1600 to 2500 are
   compared with this model's integration of the same orbit solution.
3. **ESA.** The ESA NEO Coordination Centre fits its own orbit to the observations. The JPL
   orbit is compared with it, and ESA's orbit is integrated on its own.
4. **Orbit uncertainty.** Clone orbits are drawn from the JPL covariance matrix and
   integrated alongside the nominal orbit.
5. **Thermal drift** (2025 ON9 only): the nominal orbit is integrated with a range of
   constant drifts in semi-major axis.

## 4. 2026 LO4 and Venus

### 4.1 The orbit

| | JPL solution #1 | ESA NEOCC |
|---|---|---|
| Epoch | 2026 Jun 9.0 TDB (JD 2461200.5) | 2026 Dec 26.0 (MJD 61400) |
| Semi-major axis a (au) | 0.7217597 ± 0.0000012 | 0.7217699 |
| Eccentricity e | 0.1944728 ± 0.0000025 | 0.1944425 |
| Inclination i (°) | 25.09652 ± 0.00005 | 25.09685 |
| Ascending node Ω (°) | 307.94236 ± 0.00018 | 307.93825 |
| Argument of perihelion ω (°) | 110.22967 ± 0.00032 | 110.20534 |
| Mean anomaly M (°) | 109.19708 ± 0.00052 | 70.69028 |
| Absolute magnitude H | 18.98 ± 0.29 | 19.29 |

JPL's solution is dated 21 October 2026 and uses 30 observations from 7 June to 27 September
2026 (112 days), with residuals of 0.49″ rms and condition code 4. It is JPL's first solution
for the object, which suggests the orbit became public only around that date. The object is
classed as an Atira: perihelion 0.581 au, aphelion 0.862 au, entirely inside Earth's orbit.
Its period is 223.97 days against Venus's 224.70. The closest its orbit comes to Earth's is
0.254 au. For albedos of 0.25 to 0.05 the magnitudes imply a diameter of about 0.37 to 0.95
km.

### 4.2 What it does

The semi-major axis alternates between about 0.0023 au inside Venus's and 0.0023 au outside
it, and φ swings through almost a full turn. One complete cycle takes about 330 years (332
at present). The turning points near the present are:

| Year | Angle ahead of Venus | What happens |
|---|---|---|
| 1544 | 22.6° | turns back before reaching Venus |
| 1710 | 18.1° | has passed Venus; turns back just ahead of it |
| 1877 | 23.0° | turns back before reaching Venus |
| **2043** | **18.5°** | **has passed Venus; turns back just ahead of it** |
| 2209 | 23.4° | turns back before reaching Venus |
| 2374 | 18.8° | has passed Venus; turns back just ahead of it |

So the two ends of the swing differ. Coming round from ahead of Venus, on the outer lane,
the asteroid turns back about 23° short of the planet, as in an ordinary horseshoe. Coming
up from behind, on the inner lane, it does not stop short. It passes Venus's position,
changes lane about 18° ahead of the planet, and drifts back past Venus again. Each swing
covers about 355°, and the gap it never enters is only about 5° wide, centred some 21° ahead
of Venus. This is a compound horseshoe / quasi-satellite orbit (called H-RS or HS-QS in the
literature), not a pure horseshoe.

It can pass Venus's position without a close encounter because its own orbit is tilted by
25° and moderately eccentric, which keeps it well away from the planet itself.

### 4.3 Now: a quasi-satellite phase

At the element epoch the asteroid is 4° behind Venus in mean longitude and closing. φ stays
within 25° of zero from 2014 to 2071, and during those 57 years the asteroid loops in
Venus's neighbourhood once per orbit, as a quasi-satellite does. Daily sampling of the
nominal orbit gives three groups of approaches in this phase:

| Years | Passes within 0.12 au of Venus | Closest |
|---|---|---|
| 2013.5 – 2022.0 | 15 | 0.075 au in 2017 |
| 2035.3 – 2050.0 | 25 | 0.042 au (6.3 million km) in September 2042 |
| 2063.3 – 2072.5 | 16 | 0.067 au in 2068 |

The far turns bring it a little closer still. When the asteroid stops some 23° ahead of
Venus, the near end of its loop reaches back to the planet: 0.035 au in 1877 and 0.034 au
(5.1 million km) in 2209, the closest pass in daily sampling of the seven centuries from
1700 to 2400, which span two full cycles.

Venus's sphere of gravitational control (Hill radius) is about 0.0067 au, so even those
passes are five times too far out for capture. The asteroid is never a moon. Over the same
seven centuries it comes no closer to Earth than 0.254 au (38 million km).

### 4.4 How long

In the nominal orbit the present run of regular lane switches began about 3,550 years ago
and ends about 17,500 years from now, some 21,000 years in all. Before that the asteroid
circulated slowly just outside the resonance, with at least one earlier episode of libration
between roughly 10,100 and 7,100 years ago. Afterwards it circulates just inside. It is
therefore a long-lived but temporary companion that hugs the edge of the resonance.

Over the ±20,000 years of the nominal run the eccentricity ranges from 0.14 to 0.31 and the
inclination from 28° to 21°, one rising as the other falls. Even at e = 0.31 the aphelion is
0.95 au and the perihelion 0.50 au, so in that time the orbit reaches neither Earth's nor
Mercury's. That is probably why its motion is so regular.

### 4.5 Robustness

| Check | Result |
|---|---|
| JPL epoch → JPL covariance epoch (47 days) | all six elements within 0.0001σ; perihelion direction to 4×10⁻⁸ degrees. Without relativity the error is 300 times larger. |
| JPL Horizons, 37 dates, 1600–2500 | largest difference 0.00003° in mean longitude, 2×10⁻⁹ au in a, 2×10⁻⁸ in e. Horizons' own elements show the lane switches near 1710, 1877, 2043, 2209 and 2374. |
| ESA orbit vs JPL orbit carried 200 days to ESA's epoch | all six elements agree within 0.3 of ESA's 1σ |
| ESA orbit integrated ±6,000 years | same turning points to within a year; libration begins 3,552 years ago; beside Venus 2014–2071 |
| 48 clones + nominal, ±25,000 years | libration begins 3,549–3,552 years ago and ends 17,185–17,901 years ahead in every orbit |

The 49 orbits stay within about 1° of one another in φ for 2,000 years either way, and
within about 13° at +10,000 years. By the standards of Venus co-orbitals, which the
literature describes as strongly chaotic with predictability horizons of order 150 years,
this is an unusually predictable object.

### 4.6 Why it stands out

- **Low eccentricity.** A 2025 study noted that only one of the 20 Venus co-orbitals then
  known had an eccentricity even slightly under 0.38, the value at which a Venus co-orbital's
  aphelion reaches Earth's perihelion distance. It argued that this is an observational bias
  and that a sizeable population at lower eccentricity probably exists unseen, because such
  objects never come near Earth. (A 21st was added in early 2026; it is notable for how
  closely its orbit approaches Earth's.) With e = 0.19 and an aphelion of 0.86 au, 2026 LO4
  is the kind of object that study predicted.
- **Size.** At H ≈ 19 it is a sizeable object, several hundred metres across.
- **Regularity.** More than sixty consecutive cycles in the nominal orbit, with clones that
  stay in step for millennia.

It is not the first Venus quasi-satellite (Zoozve was identified two decades ago and the
catalogues list a few others), and it is not a pure horseshoe, a type for which no Venus
example was known as of early 2026.

## 5. 2025 ON9 and Mars

### 5.1 The orbit

| | JPL solution #10 | ESA NEOCC |
|---|---|---|
| Epoch | 2026 Jun 9.0 TDB (JD 2461200.5) | same |
| a (au) | 1.5216379 ± 0.0000074 | 1.5216383 |
| e | 0.2139358 ± 0.0000026 | 0.2139359 |
| i (°) | 3.99978 ± 0.00004 | 3.99979 |
| Ω (°) | 141.98523 ± 0.00002 | 141.98524 |
| ω (°) | 198.15378 ± 0.00013 | 198.15373 |
| M (°) | 144.95091 ± 0.00105 | 144.95088 |
| H | 23.41 ± 0.26 | 23.39 |

JPL's solution (27 August 2026) uses 62 observations from 21 June to 9 November 2025 (141
days), 0.41″ rms, condition code 4. The discovery was announced in MPEC 2025-O167 on 31 July
2025. It is an Amor: perihelion 1.196 au, aphelion 1.847 au, period 685.6 days against
Mars's 687.0. Its orbit comes within 0.187 au of Earth's. The implied diameter is about 55
to 125 m.

### 5.2 What it does

A regular, nearly symmetric horseshoe. The semi-major axis alternates between about 0.0027
au inside and outside Mars's, and φ swings between 18° ahead of Mars and 18° behind it, the
long way round. One full cycle takes about 1,750 years.

| Year | Position | Lane change |
|---|---|---|
| 47 | 18.6° ahead of Mars | outer → inner |
| 903 | 18.0° behind | inner → outer |
| **1784** | **18.3° ahead** | **outer → inner** |
| 2656 | 18.4° behind | inner → outer |
| 3532 | 18.1° ahead | outer → inner |
| 4404 | 19.0° behind | inner → outer |

Today it is about 110° ahead of Mars on the inner lane and still pulling away. It will pass the
far side of the Sun from Mars around 2220 and approach Mars from behind in the 2600s.

At the turns of 1783 and 2657 it comes no closer to Mars than 0.156 au (23 million km),
more than 20 times Mars's Hill radius. Between 1690 and 2760 it comes no closer to Earth than
0.188 au (28 million km).

### 5.3 How long

| | Result |
|---|---|
| Nominal orbit, ±9,000 years | pure horseshoe throughout; e between 0.205 and 0.219, i between 2.6° and 4.2° |
| 32 clones + nominal, ±30,000 years | every orbit is a pure horseshoe from at least 8,100 years ago to at least 10,600 years ahead |
| | 25 of the 33 keep it for the entire 60,000 years |
| | the other 8 leave between 10,700 and 26,900 years ahead (the nominal orbit at 26,400); one entered 8,100 years ago |

The orbits stay within about 3° of one another for 2,000 years either way and spread to tens
of degrees by 10,000 years, so statements beyond that are statistical.

### 5.4 Checks against JPL and ESA

| Check | Result |
|---|---|
| JPL epoch → JPL covariance epoch (309 days back) | all six elements within 0.001σ |
| JPL Horizons, 37 dates, 1600–2500 | largest difference 0.0007° in mean longitude, 3×10⁻⁸ au in a |
| Horizons' own elements | 18.8° ahead of Mars and outside its orbit in 1775; 19.7° ahead and inside it in 1800: the 1784 turn, seen directly |
| ESA orbit vs JPL orbit (same epoch) | all six elements within 0.4 of ESA's 1σ |
| ESA orbit integrated ±9,000 years | same turning points to within a year |

### 5.5 Thermal drift

Sunlight re-radiated unevenly pushes small asteroids slowly inward or outward (the Yarkovsky
effect). For a body of this size near Mars a drift of order 0.001 au per million years is
typical, and neither its size nor its direction is known for 2025 ON9. The nominal orbit was
integrated ±30,000 years with constant drifts of 0, ±0.001, ±0.003 and ±0.01 au/Myr.

| Drift (au/Myr) | Pure horseshoe until |
|---|---|
| −0.01, −0.003, −0.001 | end of run (+30,000 years) |
| 0 | +26,400 years |
| +0.001 | +15,400 years |
| +0.003 | +13,900 years |
| +0.01 | +13,800 years |

In every case the horseshoe was in place for the whole 30,000 years back and survives at
least 13,800 years ahead. An outward drift shortens its future; an inward one lengthens it.

### 5.6 Context

The catalogue of Pan & Gallardo, built from orbit data of July 2024, lists five horseshoe
companions of Mars (2007 UR2, 2017 QW35, 2019 KF1, 2022 OG2, 2023 QS3) among 48 Mars
co-orbitals of all kinds. 2025 ON9 was announced a year after that data was taken. It would
be a sixth, and it resembles 2022 OG2 in eccentricity and semi-major axis. It is a routine
addition to a known class, not a first of its kind.

## 6. What the screening recovered

The screening is only as trustworthy as its record on objects whose behaviour is already
known, so that record matters more than the two new entries. Full tables are in
`results/screening/`.

| Zone | Objects | Quasi-satellite | Trojan, L4 / L5 | Horseshoe | Compound | Not in resonance now |
|---|---|---|---|---|---|---|
| Venus | 20 | 3 | 0 / 2 | 0 | 2 | 13 |
| Earth | 136 | 9 | 2 / 1 | 10 | 19 | 95 |
| Mars | 134 | 2 | 6 / 16 | 6 | 2 | 102 |

Runs covered ±3,000 years for Venus, ±1,500 for Earth and ±6,000 for Mars. A label
describes the swing the object is in at present, however short-lived that swing is.

Well-known companions, and what the screen called them:

| Object | Published description | Label here |
|---|---|---|
| 524522 Zoozve | quasi-satellite of Venus | quasi-satellite |
| 469219 Kamoʻoalewa, 164207 Cardea, 2023 FW13, 2025 PN7 | quasi-satellites of Earth | quasi-satellite |
| 706765 (2010 TK7), 614689 (2020 XL5) | Earth Trojans at L4 | L4 |
| 2002 AA29, 419624 (2010 SO16) | Earth horseshoes | horseshoe |
| 3753 Cruithne | Earth horseshoe of the compound kind | compound |
| 5261 Eureka, 101429 (1998 VF31), 311999 (2007 NS2) | Mars Trojans at L5 | L5 |
| 121514 (1999 UJ7) | Mars Trojan at L4 | L4 |
| 391595 (2007 UR2), 2017 QW35, 2019 KF1, 2022 OG2 | Mars horseshoes (Pan & Gallardo) | horseshoe |
| 322756 (2001 CK32) | Venus horseshoe–quasi-satellite compound | **not in resonance** (see below) |

Three things qualify this record.

- **One miss.** 2001 CK32 is a recognised Venus co-orbital. In this run its φ circulates up
  to the present and first reverses about seven years after the epoch; from then on it
  swings through about 358°, turning 42° and 46° behind Venus. That is the same kind of
  state as 2026 LO4's, mirrored to the trailing side. Because the present falls just before
  the first reversal, the classifier reports it as not in resonance. The screen is blunt at
  such edges, and its 360° cut-off is a convenience, not a physical boundary: 2026 LO4's own
  swing (about 355°) lies just inside it.
- **Absent objects.** The cut on orbit quality removes some known companions before the
  screen ever sees them. 2013 ND15 (a Venus Trojan), 2025 TV10 (the 21st Venus co-orbital)
  and 2023 QS3 (one of Pan & Gallardo's five Mars horseshoes) all have condition code 7 in
  the database and are not in the lists.
- **Loose labels.** "Quasi-satellite" here means any swing centred within 35° of the
  planet, including very wide ones (2017 XG62 swings 150° either side of Mars) that are
  not quasi-satellites in the usual sense. Some labelled states last only decades: 2019 EH is
  "L4" for 47 years, on an orbit fitted to 58 days of observations.

Nothing new is claimed for Earth. The Earth zone contains several recently designated
objects in or near resonance, but those checked were already in the literature (for example
2025 PN7, 2025 FP4 and 2025 SC as quasi-satellites), and the rest were not followed up.

One Earth object is worth a remark. 2005 QQ87 is at present in a wide swing confined to the
trailing side of Earth (φ between about −36° and −169°), which is Trojan-like motion about
L5, and no L5 Earth Trojan is recognised. But it is highly inclined, the episode is
temporary, and the Pan & Gallardo catalogue already lists the object as a co-orbital that
jumps between states, so this is not presented as a finding.

## 7. Was this already known?

What was checked on 10 November 2026:

- **Web searches** for each designation together with "co-orbital", "quasi-satellite",
  "horseshoe" and the planet's name returned nothing relevant. 2026 LO4 had no Wikipedia
  article and Wikipedia's list of Atira asteroids ended before it; Wikipedia's article on
  Mars trojans did not mention 2025 ON9.
- **Catalogues.** Pan & Gallardo (data of July 2024) and the Celestial Mechanics paper of
  February 2026 (data of 5 May 2025, Mars list taken from Pan & Gallardo) both predate 2025
  ON9, announced on 31 July 2025. Both predate 2026 LO4 by more than a year.
- **The latest Venus paper found** (A&A, February 2026) counts 21 Venus co-orbitals, the
  21st being 2025 TV10.
- **Timing.** JPL's first orbit solution for 2026 LO4 is dated 21 October 2026, three weeks
  before this analysis.

What this does not cover: papers in preparation or in review, conference abstracts, mailing
lists such as the Minor Planet Mailing List, social media, and anything a search engine had
not indexed. Specialists watch for exactly these objects, and a new Atira with Venus's
orbital period is conspicuous. It is quite possible that 2026 LO4's state has been noticed
and simply not yet published.

## 8. Limitations

- **No expert review.** The work was done by an AI model in a single session.
- **Short observational arcs.** 112 and 141 days. The formal uncertainties are small and the
  clones agree, but formal uncertainties assume all observations are correctly linked and
  correctly weighted. New observations will move both orbits by amounts that could change
  the long-range numbers, though probably not the present state.
- **Long-range dates are model-dependent.** Start and end dates of the resonance thousands
  of years away depend on effects left out of the model, above all thermal forces. Their
  influence was tested for 2025 ON9 only. 2026 LO4 sits close to the edge of the resonance,
  so its dates may be the more sensitive of the two.
- **Transcribed data.** Input values were relayed through a web tool and typed into files.
  The two orbits are verified by the checks; the zone lists are not verified row by row.
- **Screening is coarse and incomplete.** Objects with poorly determined orbits (condition
  code above 5) were excluded, classification labels are approximate, and only the Venus,
  Earth and Mars zones were examined.
- **DE421 rather than DE441.** Planetary initial conditions come from an older ephemeris than
  JPL uses now. The Horizons comparison shows the difference is negligible here.

## 9. What would settle it

1. An independent integration by someone who does this professionally, using their own
   software and the current orbit solutions.
2. More astrometry. 2026 LO4 never strays more than about 60° from the Sun in the sky, so
   it needs twilight observations; archival images from before June 2026 would lengthen the
   arc most.
3. A check with specialists for prior identification. The authors of the Venus co-orbital
   papers cited below, or the Minor Planet Mailing List, would know.

## 10. References and data sources

- JPL Small-Body Database, <https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html>, and its API,
  <https://ssd-api.jpl.nasa.gov/doc/sbdb.html>
- JPL Horizons, <https://ssd.jpl.nasa.gov/horizons/>
- ESA NEO Coordination Centre, <https://neo.ssa.esa.int>
- MPEC 2025-O167: 2025 ON9, <https://www.minorplanetcenter.net/mpec/K25/K25OG7.html>
- N. Pan & T. Gallardo, *An attempt to build a dynamical catalog of Solar System
  co-orbitals*, [arXiv:2410.20015](https://arxiv.org/abs/2410.20015)
- *The invisible threat: assessing the collisional hazard posed by the undiscovered Venus
  co-orbital asteroids* (2025), [arXiv:2505.15968](https://arxiv.org/abs/2505.15968)
- *The dynamically hazardous asteroid 2025 TV10: a new co-orbital asteroid of Venus*,
  A&A 706, L9 (2026), [arXiv:2601.21554](https://arxiv.org/abs/2601.21554)
- V. Carruba, G. Caritá, S. Aljbaae, R. A. N. Araujo, R. C. Domingos et al., *Co-orbital
  asteroids of terrestrial planets affected by the von Zeipel–Lidov–Kozai mechanism*,
  Celest. Mech. Dyn. Astron. 138, 8 (2026),
  [doi:10.1007/s10569-026-10275-5](https://link.springer.com/article/10.1007/s10569-026-10275-5)
- *New Transient Co-orbital Asteroids of Venus*, Res. Notes AAS (2024),
  [doi:10.3847/2515-5172/ad7261](https://iopscience.iop.org/article/10.3847/2515-5172/ad7261)
- H. Rein & S.-F. Liu, *REBOUND*, A&A 537, A128 (2012); H. Rein & D. S. Spiegel, *IAS15*,
  MNRAS 446, 1424 (2015); D. Tamayo et al., *REBOUNDx*, MNRAS 491, 2885 (2020)
