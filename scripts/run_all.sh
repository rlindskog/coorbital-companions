#!/usr/bin/env sh
# Reproduce the whole analysis.
#
#   sh scripts/run_all.sh            full run, as in the write-up: about 2 hours on 2 CPU cores
#   QUICK=1 sh scripts/run_all.sh    short spans, about 7 minutes; writes to quick/ and leaves
#                                    the committed results/ and figures/ untouched
set -e
cd "$(dirname "$0")/.."
PY=${PYTHON:-python3}

if [ "${QUICK:-0}" = "1" ]; then
  export COORB_RESULTS_DIR=quick/results COORB_FIGURES_DIR=quick/figures COORB_RUNS_DIR=quick/runs
  SCREEN_V="--years 600"; SCREEN_E="--years 300"; SCREEN_M="--years 1500"
  NOM_LO4="--years 1000"; NOM_ON9="--years 2800"
  ESA="--years 1000"; CLONES="--n 8 --years 1000"; DRIFT="--years 1000"
else
  SCREEN_V=""; SCREEN_E=""; SCREEN_M=""; NOM_LO4=""; NOM_ON9=""; ESA=""; CLONES=""; DRIFT=""
fi

# 1. discovery step: screen each planet's co-orbital zone
$PY scripts/screen_zone.py venus $SCREEN_V
$PY scripts/screen_zone.py earth $SCREEN_E
$PY scripts/screen_zone.py mars  $SCREEN_M

# 2. checks of the model against JPL and ESA, and the nominal orbits
for OBJ in 2026_LO4 2025_ON9; do
  $PY scripts/check_propagation.py $OBJ
  $PY scripts/compare_horizons.py  $OBJ
  $PY scripts/close_approaches.py  $OBJ
  $PY scripts/compare_esa.py       $OBJ $ESA
done
$PY scripts/nominal.py 2026_LO4 $NOM_LO4
$PY scripts/nominal.py 2025_ON9 $NOM_ON9

# 3. robustness: orbit uncertainty (clones) and thermal drift
$PY scripts/clones.py 2026_LO4 $CLONES
$PY scripts/clones.py 2025_ON9 $CLONES
$PY scripts/yarkovsky.py 2025_ON9 $DRIFT

# 4. figures
$PY scripts/make_figures.py
