#!/bin/bash
# ATLAS-vs-CMS comparison of a production at scale: one array task per seed (all checks + a slim summary), then one
# merge job (check tally + plots). From the repo root, with the site options of slurm/submit.sbatch:
#
#   BASE=<FINAL_BASE of the production, holding seed*/>  OUT=<directory for summaries/, logs and plots>
#   N=$(ls -d $BASE/seed* | wc -l)
#   jid=$(BASE=$BASE OUT=$OUT sbatch --parsable --array=0-$((N-1)) slurm/compare.sbatch)
#   BASE=$BASE OUT=$OUT sbatch --dependency=afterany:$jid --mem=32G --time=04:00:00 slurm/compare.sbatch   # no --array: merge
#SBATCH --job-name=mdcmp
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
set -eo pipefail
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
[ -f "${REPO_DIR}/scripts/compare.sh" ] || { echo "submit from the repo root, or set REPO_DIR"; exit 1; }
: "${BASE:?directory holding seed*/}" "${OUT:?output directory}"
mkdir -p "${OUT}/summaries"
# the regenerated HepMC text goes to TMPDIR: not the RAM-backed /tmp of a diskless node
export TMPDIR="${SCRATCH_BASE:-${PSCRATCH:-${TMPDIR:-/tmp}}}/mdcmp_${USER}_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-m}"
mkdir -p "${TMPDIR}"; trap 'rm -rf "${TMPDIR}"' EXIT
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  DIRS=($(ls -d "${BASE}"/seed*)); d="${DIRS[${SLURM_ARRAY_TASK_ID}]}"; s="$(basename "${d}")"
  "${REPO_DIR}/scripts/compare.sh" --no-plots --summary "${OUT}/summaries/${s}.pkl" "${d}" > "${OUT}/summaries/${s}.log" 2>&1 || rc=$?
  tail -1 "${OUT}/summaries/${s}.log"; exit ${rc:-0}
else
  "${REPO_DIR}/scripts/compare.sh" --merge --out "${OUT}" --title "${TITLE:-$(basename "${BASE}")}" "${OUT}"/summaries/*.pkl 2>&1 | tee "${OUT}/compare.log"
fi
