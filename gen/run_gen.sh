#!/bin/bash
# Step 1: generate the common Pythia8 event sample (HepMC3 + HepMC2 ASCII).
#
#   gen/run_gen.sh <process> <seed> <nevents> <outdir>
#   e.g. gen/run_gen.sh zee_13p6TeV 12345 1000 output/zee_13p6TeV_seed12345
#
# Output: <outdir>/events.hepmc3, <outdir>/events.hepmc2, <outdir>/events.json, <outdir>/gen.log
# Reproducible: same card + seed + Pythia version  =>  byte-identical output (verified by gen/test_reproducible.sh).
set -eo pipefail
PROC="${1:?process name, e.g. zee_13p6TeV (see config/pythia/)}"
SEED="${2:?random seed (int)}"
NEV="${3:?number of events}"
OUT="${4:?output directory}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO}/scripts/common.sh"
CARD="${REPO}/config/pythia/${PROC}.cmnd"
[ -f "${CARD}" ] || { echo "no such card: ${CARD}"; exit 1; }
# The LCG view is an EL9 build: on any other host OS (e.g. SLES on Perlmutter) re-run this script inside the EL9 image.
if [ -z "${GEN_IN_CONTAINER:-}" ] && ! grep -qs 'platform:el9' /etc/os-release; then
  source "${REPO}/config/versions.env"
  mkdirp "${OUT}"; OUT="$(cd "${OUT}" && pwd)"
  GEN_IN_CONTAINER=1 exec apptainer exec --no-home $(apptainer_binds) -B "${REPO}" -B "${OUT}" "${ATLAS_CONTAINER}" \
    bash "${REPO}/gen/run_gen.sh" "${PROC}" "${SEED}" "${NEV}" "${OUT}"
fi
[ -x "${REPO}/gen/build/gen_hepmc" ] || bash "${REPO}/gen/build.sh"
source "${REPO}/gen/env_lcg.sh"
mkdirp "${OUT}"; OUT="$(cd "${OUT}" && pwd)"
cp "${CARD}" "${OUT}/pythia.cmnd"
echo "[gen] ${PROC} seed=${SEED} nev=${NEV} -> ${OUT}"
( cd "${OUT}" && "${REPO}/gen/build/gen_hepmc" pythia.cmnd "${SEED}" "${NEV}" events ) > "${OUT}/gen.log" 2>&1
tail -1 "${OUT}/gen.log"
sha256sum "${OUT}/events.hepmc3" "${OUT}/events.hepmc2" | tee "${OUT}/SHA256SUMS"
