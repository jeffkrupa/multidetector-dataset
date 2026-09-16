#!/bin/bash
# End-to-end driver: Pythia8 -> {ATLAS, CMS} full simulation + reconstruction of the SAME events.
#
#   ./run_all.sh <process> <seed> <nevents> [outdir]
#   e.g. ./run_all.sh zee_13p6TeV 12345 100
#
# Layout of <outdir> (default: output/<process>_seed<seed>_n<nevents>):
#   gen/    events.hepmc3, events.hepmc2, events.json, pythia.cmnd, SHA256SUMS
#   atlas/  EVNT -> HITS -> AOD -> DAOD_PHYSLITE (+ logs)
#   cms/    GEN-SIM -> DIGI-RAW-HLT -> AOD/MiniAOD -> NanoAOD (+ logs)
#
# Each step is a standalone script that can be re-run on its own; the experiment steps
# run in parallel (they are independent given the HepMC file). Set STEPS to restrict, e.g.
#   STEPS="gen cms" ./run_all.sh ...
set -eo pipefail
PROC="${1:?process name, see config/pythia/}"; SEED="${2:?seed}"; NEV="${3:?nevents}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${REPO}/scripts/common.sh"
OUT="${4:-${REPO}/output/${PROC}_seed${SEED}_n${NEV}}"
STEPS="${STEPS:-gen atlas cms}"
mkdirp "${OUT}"; OUT="$(cd "${OUT}" && pwd)"
echo "[run_all] process=${PROC} seed=${SEED} nevents=${NEV} steps='${STEPS}'"
echo "[run_all] output: ${OUT}"

if [[ " ${STEPS} " == *" gen "* ]]; then
  "${REPO}/gen/run_gen.sh" "${PROC}" "${SEED}" "${NEV}" "${OUT}/gen"
fi

pids=()
if [[ " ${STEPS} " == *" atlas "* ]]; then
  ( "${REPO}/atlas/run_atlas.sh" "${OUT}/gen/events.hepmc3" "${OUT}/atlas" "${SEED}" "${NEV}" > "${OUT}/atlas.log" 2>&1 \
      && echo "[run_all] ATLAS done" || { echo "[run_all] ATLAS FAILED, see ${OUT}/atlas.log"; exit 1; } ) & pids+=($!)
fi
if [[ " ${STEPS} " == *" cms "* ]]; then
  ( "${REPO}/cms/run_cms.sh" "${OUT}/gen/events.hepmc2" "${OUT}/cms" "${SEED}" "${NEV}" > "${OUT}/cms.log" 2>&1 \
      && echo "[run_all] CMS done" || { echo "[run_all] CMS FAILED, see ${OUT}/cms.log"; exit 1; } ) & pids+=($!)
fi
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
[ $rc -eq 0 ] && echo "[run_all] all steps finished: ${OUT}"
exit $rc
