#!/bin/bash
# Create the CMSSW project area once (inside the CMS container) and compile the HepMC relabeler plugin.
# Output: ${CMS_WORK_AREA:-cms/build}/${CMSSW_RELEASE}
set -eo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO}/config/versions.env"; source "${REPO}/scripts/common.sh"
WORK="${CMS_WORK_AREA:-${REPO}/cms/build}"; mkdirp "${WORK}"; WORK="$(cd "${WORK}" && pwd)"
cat > "${WORK}/build_inner.sh" <<INNER
set -eo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=${SCRAM_ARCH}
cd "${WORK}"
[ -d ${CMSSW_RELEASE}/src ] || scram project CMSSW ${CMSSW_RELEASE}
cd ${CMSSW_RELEASE}/src && eval \$(scram runtime -sh)
rm -rf HepMCTools && cp -r "${REPO}/cms/HepMCTools" .
scram b -j 4 > "${WORK}/scram_build.log" 2>&1 || { tail -30 "${WORK}/scram_build.log"; exit 1; }
echo "[cms] built ${CMSSW_RELEASE} (${SCRAM_ARCH}) + HepMCTools/Relabel in ${WORK}/${CMSSW_RELEASE}"
INNER
exec "${CMS_CONTAINER_WRAPPER}" --no-home --pwd "${WORK}" $(apptainer_binds) -B "${REPO}" -- bash "${WORK}/build_inner.sh"
