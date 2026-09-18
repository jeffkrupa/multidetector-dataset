#!/bin/bash
# Slurm job array: N independent jobs, one seed each, through condor/job.sh (the wrapper is scheduler-agnostic).
#
#   condor/make_tarball.sh                                   # pack the repo state to ship to the nodes
#   sbatch --array=0-19 slurm/submit.sbatch                  # 20 jobs, seeds seed0..seed0+19
#   PROCESS=zee_13p6TeV NEVENTS=200 SEED0=2000 KEEP=min HEPMC=0 sbatch --array=0-9 --cpus-per-task=1 slurm/submit.sbatch
#
# Requirements on the nodes: x86_64 EL9, CVMFS (sft, atlas, atlas-condb, cms, unpacked), apptainer, outbound access to the
# ATLAS and CMS Frontier conditions servers. Outputs are written to $FINAL_BASE on a filesystem the nodes can write.
#SBATCH --job-name=multidet
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --tmp=8G
#SBATCH --output=%x_%A_%a.out
set -eo pipefail
REPO_DIR="${REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PROCESS="${PROCESS:-ttbar_13p6TeV}"; NEVENTS="${NEVENTS:-100}"; SEED0="${SEED0:-1000}"; KEEP="${KEEP:-min}"; HEPMC="${HEPMC:-0}"
SEED=$((SEED0 + ${SLURM_ARRAY_TASK_ID:-0}))
FINAL_BASE="${FINAL_BASE:-${REPO_DIR}/output/slurm/${PROCESS}_n${NEVENTS}}"
export REPO_TARBALL="${REPO_DIR}/condor/repo.tar.gz"     # shipped state; the wrapper unpacks and builds it on the node
export SLURM_TMPDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp/${USER}_${SLURM_JOB_ID:-$$}}}"
exec bash "${REPO_DIR}/condor/job.sh" "${PROCESS}" "${SEED}" "${NEVENTS}" "${FINAL_BASE}/seed${SEED}" "${KEEP}" "${HEPMC}"
