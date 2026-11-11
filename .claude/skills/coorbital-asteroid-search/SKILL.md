---
name: coorbital-asteroid-search
description: Find and validate unreported co-orbital asteroids (Trojans, horseshoes, quasi-moons) of a planet from NASA/JPL orbit data, or check whether a given asteroid is one.
---

# Co-orbital asteroid search

Screen catalogued asteroids near a planet's orbit, simulate each one, and pick out those locked in the planet's 1:1 resonance that the literature has not yet described. This procedure found 2026 LO4 (Venus, compound horseshoe and quasi-satellite) and 2025 ON9 (Mars, horseshoe) in November 2026.

## What can be claimed

The asteroids are already discovered by sky surveys. What this finds is how one moves relative to a planet, which nobody may have worked out yet. Hold to these rules in everything you write:

- Say "apparently unreported" or "newly recognised", never "discovered an asteroid" or "confirmed new". A literature search cannot prove nobody has noticed.
- Never call the object a moon. A quasi-satellite orbits the Sun; the loop around the planet appears only from a viewpoint moving with the planet, and the phase is temporary.
- Check every counting claim ("the first", "the second", "the only") against the literature on its own. Hybrids are easy to miss: Venus already has horseshoe and quasi-satellite hybrids (2001 CK32, 2024 AF6).
- Give each number with its window, for example "closest pass to Venus between 1700 and 2400".
- Always state the observation count, arc length and condition code, and that the work is AI-assisted and unreviewed until a specialist checks it.
- A screen that finds nothing new is a normal result. Report the census and stop. Do not stretch a weak case.

## Why this works

New asteroids are designated every week. Published co-orbital catalogues lag by a year or more. An object first observed after a catalogue's data cut-off, with an orbit good enough to simulate, has usually not been classified by anyone. That gap is where to look.

## Setup

This skill ships inside the repository it drives. Work from the repository root. If you are somewhere else, clone it first:

```sh
git clone https://github.com/rlindskog/coorbital-companions
cd coorbital-companions
```

Then install and run the fast checks:

```sh
pip install -r requirements.txt     # numpy, rebound, reboundx, jplephem, de421, matplotlib
python -m pytest -q                 # six fast checks
```

Read `README.md` and section 3 of `WRITEUP.md` first. If a script is missing, rebuild from the Method section below.

## Step 1: get the candidate list

A planet's co-orbital zone is its semi-major axis plus or minus `(m_planet / 3 M_sun)^(1/3)` times itself.

| Planet | a range (au) | Screen span each way | Cadence | Smoothing |
|---|---|---|---|---|
| Venus | 0.7165 to 0.7301 | 3,000 yr | 0.1 yr | 2 yr |
| Earth | 0.990 to 1.010 | 1,500 yr | 0.2 yr | 3 yr |
| Mars | 1.5164 to 1.5310 | 6,000 yr | 0.5 yr | 5 yr |

Only these three zones are tested. For another planet, compute the range from the formula and add an entry to `ZONES` in `scripts/screen_zone.py` and `scripts/fetch_data.py`: span of at least ten libration periods, cadence about a quarter of the planet's period or less, smoothing about three planet periods.

Query the JPL Small-Body Database for every object in the range with condition code 5 or better:

```
https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=pdes,epoch,a,e,i,om,w,ma,condition_code,data_arc,first_obs&sb-cdata={"AND":["a|RG|LO|HI","condition_code|LE|5"]}&sort=a&full-prec=1
```

- Try `python scripts/fetch_data.py zone venus` first. It writes to `data/fresh/`. It was written without access to the live API, so check its output against the raw response it saves.
- If the shell cannot reach the host, fetch the same URL with the web fetch tool and write the rows to a CSV by hand. Ask for the rows verbatim at full precision, split the range if the list is long, and check the row count against the response's `count`. Do not try other routes around a blocked host.
- The screen reads `data/zones/ZONE.csv` with this header: `designation,epoch_jd_tdb,a_au,e,i_deg,node_deg,peri_deg,M_deg,condition_code,data_arc_days,first_obs`. Keep the old snapshot under a dated name before replacing it, and record the retrieval date. Orbit solutions change as observations arrive.

## Step 2: screen the zone

```sh
python scripts/screen_zone.py venus          # about 3 min on 2 cores; earth 4, mars 12
python scripts/screen_zone.py venus --classify-only   # redo labels from the saved run
```

This writes `results/screening/ZONE.csv` with a label per object: QS, L4, L5, HS, compound, edge, passing, unstable.

Before trusting it, confirm it recovers known companions: Zoozve as QS (Venus); Kamoʻoalewa as QS, 2010 TK7 and 2020 XL5 as L4, Cruithne as compound (Earth); Eureka as L5, 1999 UJ7 as L4 (Mars). If these fail, fix the run before looking for anything new.

Known weaknesses of the screen:

- It labels only the swing the object is in now. An object whose first reversal falls just after the present reads as "passing" or "edge" (this happened to 2001 CK32). Look at "edge" rows too.
- The condition-code cut removes some known companions before the screen sees them.
- "QS" means any swing centred within 35 degrees of the planet, including very wide ones.
- Labels on short arcs can describe states lasting only decades.

## Step 3: pick candidates

Keep objects that are librating now, stay librating for several cycles, and have a usable orbit. Rank first those whose `first_obs` is later than the newest catalogue's data cut-off.

Then check whether each is already known. Record what you checked and what you could not:

- Pan and Gallardo's catalogue of co-orbitals (arXiv:2410.20015, data of July 2024).
- The newest papers on that planet's co-orbitals. For Venus: arXiv:2505.15968, arXiv:2601.21554. Papers by C. and R. de la Fuente Marcos often report new ones first.
- Web searches for the designation with "co-orbital", "quasi-satellite", "horseshoe", "Trojan" and the planet's name.
- Wikipedia's pages on the object, on quasi-satellites, horseshoe orbits and the planet's trojans.

Not covered by any of this: papers in review, conference abstracts, mailing lists. Say so.

## Step 4: scrutinise each survivor

Create `data/objects/NAME.json` using `data/objects/2026_LO4.json` as the template. It needs:

- `jpl`: elements at the standard epoch plus the covariance block, from `https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=DESIGNATION&cov=mat&phys-par=1&full-prec=1`. Covariance labels must be in the order e, q, tp, node, peri, i.
- `esa`: the independent orbit from `https://neo.ssa.esa.int/PSDB-portlet/download?file=DESIGNATION_NO_SPACE.ke1`.
- `horizons`: osculating elements every 25 years from 1600 to 2500, saved to `data/horizons/NAME.csv` (the URL pattern is in the template).
- `planet` and an `analysis` block with spans, cadences, clone count and the close-approach window.

`python scripts/fetch_data.py object "2026 LO4"` attempts all three downloads.

Run the checks in this order. Stop and investigate if one fails.

| Script | Question | Pass looks like |
|---|---|---|
| `check_propagation.py NAME` | Does the model carry JPL's elements at one epoch onto JPL's elements at the covariance epoch? | Every element within about 0.01 sigma. A larger gap usually means a transcription, epoch or frame error. |
| `compare_horizons.py NAME` | Does it follow JPL Horizons from 1600 to 2500? | Mean longitude within about 0.001 degrees, semi-major axis within 1e-7 au. |
| `compare_esa.py NAME` | Does ESA's own orbit behave the same? | Same state, similar turning points and timing. |
| `nominal.py NAME` | Cycle length, turning points, how long the resonance lasts | Many consecutive cycles around the present. |
| `close_approaches.py NAME` | Real distances to the planet and to Earth | Use a window of several centuries each way. A narrow window gave a wrong "closest pass" once. |
| `clones.py NAME` | Do orbits drawn from the covariance agree? | All clones librating now. Report the minimum, median and maximum of the start and end dates. |
| `yarkovsky.py NAME` | Does thermal drift break the state? | Run for small objects, roughly H above 21. |

Long runs take 20 to 30 minutes each. Run them in the background with a log file. Each script accepts `--years` and `--n` for a quick first pass.

Two extra checks the repository has no script for:

- **Secular behaviour.** Compute the argument of perihelion measured in the planet's orbital plane, and eccentricity against inclination, over the long run. A wide swing of that angle around 0 or 180 degrees indicates von Zeipel, Lidov and Kozai type evolution, which specialists will ask about.
- **Size.** Diameter in km is `1329 / sqrt(albedo) * 10^(-H/5)`. Quote the range for albedo 0.05 to 0.25.

## Method, in brief

Use this to rebuild if a script is unavailable.

- **Model.** Sun and eight planets as point masses, asteroids as massless particles, REBOUND's IAS15 integrator, units of au and days. Planet positions and masses at the asteroid's epoch come from JPL DE421 through `jplephem`, rotated to ecliptic J2000. Asteroid elements are heliocentric ecliptic J2000. Screening is Newtonian with Earth and Moon merged. Scrutiny adds the Sun's relativistic term (REBOUNDx `gr`) and, for comparisons with JPL, a separate Moon. Integrate backward and forward from the epoch.
- **Resonant angle.** phi is the asteroid's mean longitude minus the planet's. Unwrap it, smooth it over a few planet periods, and find its turning points.
- **Labels.** phi circulating means no resonance. A swing around 0 is a quasi-satellite, around plus or minus 60 a Trojan, around 180 without reaching 0 a horseshoe. One swing that encloses both 0 and 180 is compound. The librating interval is the stretch around the present in which every swing between turning points is under 360 degrees.
- **Clones.** Draw from a multivariate normal in JPL's covariance parameters at the covariance epoch and convert to elements.

## Reporting

Deliver a short table (orbit, type, cycle length, current phase and its dates, how long the resonance lasts across clones, size), one figure (phi against time, and the path in the frame rotating with the planet), and a write-up with a limitations section. `WRITEUP.md` is the model to follow.

State plainly in the limitations: short arcs, long-range dates that depend on effects left out of the model, any data transcribed by hand, and the coarse screen.

## After a finding

Suggest the user write to a specialist, ideally an author of the newest paper on that planet's co-orbitals. The email should give the elements, the method, the checks, links to code, the AI disclosure, and two questions: does this look right, and is it already known. If it holds up, Research Notes of the AAS is the usual venue for a short report. Check its current word limit and its policy on AI use before drafting; an AI cannot be an author.
