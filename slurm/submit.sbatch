#!/bin/bash
# Slurm job array: N independent jobs, one seed each, through condor/job.sh (the wrapper is scheduler-agnostic).
#
#   condor/make_tarball.sh                                   # pack the repo state to ship to the nodes
#   sbatch --array=0-19 slurm/submit.sbatch                  # 20 jobs, seeds seed0..seed0+19
#   PROCESS=zee_13p6TeV NEVENTS=200 SEED0=2000 KEEP=min HEPMC=0 sbatch --array=0-9 --cpus-per-task=1 slurm/submit.sbatch
#
# Requirements on the nodes: x86_64, CVMFS (sft, atlas, atlas-condb, cms, unpacked, oasis if no apptainer), outbound access to the
# ATLAS and CMS Frontier conditions servers. Outputs are written to $FINAL_BASE on a filesystem the nodes can write.
#SBATCH --job-name=multidet
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.out
# Site options go on the command line; a node-local disk request (--tmp=8G) where the nodes have one. NERSC Perlmutter:
#   sbatch -A <account> -C cpu -q shared --array=0-19 slurm/submit.sbatch      # from the repo root
set -eo pipefail
# Slurm runs a spooled copy of this script, so its own path says nothing: take the repo from the submit directory.
REPO_DIR="${REPO_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
[ -f "${REPO_DIR}/run_all.sh" ] || { echo "submit from the repo root, or set REPO_DIR"; exit 1; }
PROCESS="${PROCESS:-ttbar_13p6TeV}"; NEVENTS="${NEVENTS:-100}"; SEED0="${SEED0:-1000}"; KEEP="${KEEP:-min}"; HEPMC="${HEPMC:-0}"
SEED=$((SEED0 + ${SLURM_ARRAY_TASK_ID:-0}))
FINAL_BASE="${FINAL_BASE:-${REPO_DIR}/output/slurm/${PROCESS}_n${NEVENTS}}"
export REPO_TARBALL="${REPO_DIR}/condor/repo.tar.gz"     # shipped state; the wrapper unpacks and builds it on the node
# scratch: always a directory of this task's own (tasks sharing a node must not share repo/ and keep/), under the site's
# per-job dir if it sets one, else $PSCRATCH (Perlmutter: no local disk, and TMPDIR=/tmp is RAM charged to the job), else TMPDIR.
SCRATCH_BASE="${SCRATCH_BASE:-${SLURM_TMPDIR:-${PSCRATCH:-${TMPDIR:-/tmp}}}}"
export SLURM_TMPDIR="${SCRATCH_BASE}/multidet_${USER}_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"
rc=0; bash "${REPO_DIR}/condor/job.sh" "${PROCESS}" "${SEED}" "${NEVENTS}" "${FINAL_BASE}/seed${SEED}" "${KEEP}" "${HEPMC}" || rc=$?
# the task's directory is ours alone: drop what is left of it (tool temp files), unless delivery failed (rc >= 100)
[ "${rc}" -lt 100 ] && rm -rf "${SLURM_TMPDIR}"
exit ${rc}
