#!/bin/bash
# Step 2a: ATLAS full Geant4 simulation + reconstruction of an external HepMC file.
#
#   atlas/run_atlas.sh <events.hepmc3|hepmc2> <outdir> <seed> <nevents>
#
# Chain (Athena ${ATHENA_RELEASE}, Run-3 MC23 geometry/conditions, NO pile-up, no trigger):
#   HepMC --Gen_tf.py--> EVNT --Sim_tf.py (FullG4MT_QS)--> HITS --Reco_tf.py--> AOD --Derivation_tf.py--> DAOD_PHYSLITE
#   The single DAOD_PHYSLITE also carries the constituent-level content (particle-flow objects before overlap removal,
#   full truth, all tracks, topoclusters) added by atlas/python/MultiDetPHYSLITE.py (ATLAS_EXTRA_CONTENT).
# Runs inside the ATLAS AlmaLinux9 apptainer image with CVMFS; conditions come from Frontier (no grid proxy needed).
# Every step writes log.<step> and a full transform report in <outdir>.
set -eo pipefail
HEPMC="${1:?input HepMC file}"; OUT="${2:?output dir}"; SEED="${3:?seed}"; NEV="${4:?nevents}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO}/scripts/common.sh"
source "${REPO}/config/versions.env"
mkdirp "${OUT}"; OUT="$(cd "${OUT}" && pwd)"; HEPMC="$(readlink -f "${HEPMC}")"

# ---- physics/conditions configuration (edit here, not in the commands below) ----
ATLAS_GEOMETRY="${ATLAS_GEOMETRY:-ATLAS-R3S-2021-03-02-00}"
ATLAS_CONDITIONS="${ATLAS_CONDITIONS:-OFLCOND-MC23-SDR-RUN3-11-02}"
ATLAS_SIMULATOR="${ATLAS_SIMULATOR:-FullG4MT_QS}"
ATLAS_SIM_PREINCLUDE="${ATLAS_SIM_PREINCLUDE:-Campaigns.MC23eSimulationMultipleIoV}"
ATLAS_RECO_PREINCLUDE="${ATLAS_RECO_PREINCLUDE:-Campaigns.MC23eNoPileUp}"
ATLAS_ECM_GEV="${ATLAS_ECM_GEV:-13600}"
ATLAS_NTHREADS="${ATLAS_NTHREADS:-4}"   # AthenaMT for Sim_tf/Reco_tf; output is identical to 1 thread (verified), sim wall time 2.4x faster
ATLAS_DERIV_FORMATS="${ATLAS_DERIV_FORMATS:-PHYSLITE}"         # add FTAG1 for the official flavour-tagging format (post-overlap-removal FlowElements)
ATLAS_EXTRA_CONTENT="${ATLAS_EXTRA_CONTENT:-pflow truth tracks clusters}"   # constituent-level content added to DAOD_PHYSLITE (atlas/python/MultiDetPHYSLITE.py); "" = plain PHYSLITE
ATLAS_PFLOW_DUMP="${ATLAS_PFLOW_DUMP:-0}"                        # 1: also write the flat pflow.root via atlas/dump_pflow.py (same content as the pflow+truth groups above)
export ATLAS_GEOMETRY ATLAS_CONDITIONS ATLAS_SIMULATOR ATLAS_SIM_PREINCLUDE ATLAS_RECO_PREINCLUDE ATLAS_ECM_GEV ATLAS_NTHREADS ATLAS_DERIV_FORMATS ATLAS_EXTRA_CONTENT ATLAS_PFLOW_DUMP

# The inner script runs inside the container. It is written to the output dir so the exact
# commands used are preserved next to the results.
cat > "${OUT}/atlas_chain.sh" <<INNER
#!/bin/bash
set -eo pipefail
cd "${OUT}"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
set +e   # atlasLocalSetup/asetup do not survive errexit
source \${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh -q
asetup Athena,${ATHENA_RELEASE}
set -e
command -v Sim_tf.py >/dev/null || { echo "[atlas] asetup Athena,${ATHENA_RELEASE} failed"; exit 1; }
unset X509_USER_PROXY
export ATHENA_CORE_NUMBER=${ATLAS_NTHREADS}
echo "[atlas] Athena ${ATHENA_RELEASE} \$(cat \$AtlasVersion 2>/dev/null) on \$(hostname), \$(date)"

# Gen_tf.py locates the input as glob('*<name>.*ev*ts') in the CWD -> must be named *.events and live here.
ln -sf "${HEPMC}" ./input.events
rm -rf ./999002 && cp -r "${REPO}/atlas/jobconfig/999002" ./999002

step() { local name=\$1; shift; echo "[atlas] === \$name: \$(date) ==="; echo "\$@" > "cmd.\$name"; /usr/bin/time -v "\$@" > "log.\$name" 2>&1 || { echo "[atlas] \$name FAILED (see ${OUT}/log.\$name)"; tail -40 "log.\$name"; exit 1; }; grep -E "Elapsed \(wall|Maximum resident" "log.\$name" | sed "s/^/[atlas] \$name /"; }

[ -s EVNT.pool.root ] || step 1_gen  Gen_tf.py --ecmEnergy=${ATLAS_ECM_GEV} --maxEvents=${NEV} --randomSeed=${SEED} \
      --jobConfig=999002 --inputGeneratorFile=input.events --outputEVNTFile=EVNT.pool.root

[ -s HITS.pool.root ] || step 2_sim  Sim_tf.py --CA --simulator ${ATLAS_SIMULATOR} \
      --conditionsTag default:${ATLAS_CONDITIONS} --geometryVersion default:${ATLAS_GEOMETRY} \
      --preInclude "EVNTtoHITS:${ATLAS_SIM_PREINCLUDE}" --postInclude "default:PyJobTransforms.UseFrontier" \
      --inputEVNTFile EVNT.pool.root --outputHITSFile HITS.pool.root \
      --maxEvents ${NEV} --jobNumber ${SEED} --randomSeed ${SEED} --multithreaded $( [ ${ATLAS_NTHREADS} -gt 1 ] && echo True || echo False )

[ -s AOD.pool.root ] || step 3_reco Reco_tf.py --CA --inputHITSFile HITS.pool.root --outputAODFile AOD.pool.root \
      --conditionsTag default:${ATLAS_CONDITIONS} --geometryVersion default:${ATLAS_GEOMETRY} \
      --preInclude "all:${ATLAS_RECO_PREINCLUDE}" --postInclude "all:PyJobTransforms.UseFrontier" \
      --digiSeedOffset1 ${SEED} --digiSeedOffset2 ${SEED} --maxEvents ${NEV} --multithreaded $( [ ${ATLAS_NTHREADS} -gt 1 ] && echo True || echo False )

have_all=1; for f in ${ATLAS_DERIV_FORMATS}; do [ -s DAOD_\$f.out.pool.root ] || have_all=0; done
export PYTHONPATH="${REPO}/atlas/python:\${PYTHONPATH}"; export ATLAS_EXTRA_CONTENT="${ATLAS_EXTRA_CONTENT}"
EXTRA=""; [ -n "${ATLAS_EXTRA_CONTENT}" ] && EXTRA='--postInclude default:MultiDetPHYSLITE.AddMultiDetContent'
[ \$have_all = 1 ] || step 4_deriv Derivation_tf.py --CA --inputAODFile AOD.pool.root \
      --outputDAODFile out.pool.root --formats ${ATLAS_DERIV_FORMATS} --maxEvents ${NEV} \$EXTRA

[ "${ATLAS_PFLOW_DUMP}" != 1 ] || [ -s pflow.root ] || step 5_pflow python "${REPO}/atlas/dump_pflow.py" AOD.pool.root pflow.root

echo "[atlas] done: \$(date)"; ls -la *.pool.root pflow.root 2>/dev/null || true
INNER
chmod +x "${OUT}/atlas_chain.sh"

echo "[atlas] container=${ATLAS_CONTAINER} release=Athena,${ATHENA_RELEASE} input=${HEPMC} out=${OUT}"
exec apptainer exec --no-home --pwd "${OUT}" $(apptainer_binds) -B "${REPO}" -B "${OUT}" "${ATLAS_CONTAINER}" bash "${OUT}/atlas_chain.sh"
