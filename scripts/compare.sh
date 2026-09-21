#!/bin/bash
# scripts/compare_events.py with the LCG python stack (uproot, awkward, matplotlib) on any host: directly on EL9, inside
# the EL9 image elsewhere (e.g. SLES on Perlmutter). All arguments are passed through:
#
#   scripts/compare.sh output/<sample>
#   scripts/compare.sh --no-plots --summary sum/seed1000.pkl <dir>/seed1000
#   scripts/compare.sh --merge --out <plotdir> sum/*.pkl
#   scripts/compare.sh --particles --hist h/seed1000.npz <dir>/seed1000      # scripts/compare_particles.py instead
set -eo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO}/scripts/common.sh"
if [ -z "${CMP_IN_CONTAINER:-}" ] && ! grep -qs 'platform:el9' /etc/os-release; then
  source "${REPO}/config/versions.env"
  # bind the directory above every path argument (inputs exist; --out / --summary targets may not, their parent must)
  binds=$(for a in "$@"; do case "$a" in -*) continue;; esac; p="$(dirname "$(readlink -m "$a")")"; [ -d "$p" ] && echo "$p"; done | sort -u | sed 's/^/-B /')
  CMP_IN_CONTAINER=1 exec apptainer exec --no-home $(apptainer_binds) -B "${REPO}" ${binds} "${ATLAS_CONTAINER}" bash "${REPO}/scripts/compare.sh" "$@"
fi
source "${REPO}/gen/env_lcg.sh" > /dev/null 2>&1
export MPLBACKEND=Agg
SCRIPT=compare_events.py; [ "${1:-}" = "--particles" ] && { SCRIPT=compare_particles.py; shift; }   # particle-level marginals
exec python3 "${REPO}/scripts/${SCRIPT}" "$@"
