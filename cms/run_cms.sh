#!/bin/bash
# Step 2b: CMS full Geant4 simulation + reconstruction of an external HepMC2 file.
#
#   cms/run_cms.sh <events.hepmc2> <outdir> <seed> <nevents>
#
# Chain (CMSSW ${CMSSW_RELEASE}, Run-3 2024 geometry/conditions, NO pile-up, 2024 HLT menu):
#   HepMC2 --MCFileSource--> GEN-SIM --DIGI,L1,DIGI2RAW,HLT:2024v14--> GEN-SIM-RAW
#          --RAW2DIGI,L1Reco,RECO,RECOSIM,PAT--> AODSIM + MINIAODSIM --NANO(+PF candidates)--> NANOAODSIM
# Runs inside the cmssw-el9 apptainer wrapper with CVMFS, on node-local scratch (CMS_SCRATCH, default \$TMPDIR or /tmp). Conditions from the CMS Frontier/CVMFS setup done by
# cmsset_default (no proxy needed). Every step keeps its full cmsDriver config (step*_cfg.py) and log in <outdir>.
set -eo pipefail
HEPMC="${1:?input HepMC2 file}"; OUT="${2:?output dir}"; SEED="${3:?seed}"; NEV="${4:?nevents}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO}/config/versions.env"; source "${REPO}/scripts/common.sh"
mkdirp "${OUT}"; OUT="$(cd "${OUT}" && pwd)"; HEPMC="$(readlink -f "${HEPMC}")"
WORK="${CMS_WORK_AREA:-${REPO}/cms/build}"
[ -d "${WORK}/${CMSSW_RELEASE}/src/HepMCTools" ] || bash "${REPO}/cms/build.sh"
WORK="$(cd "${WORK}" && pwd)"

# ---- physics/conditions configuration (edit here, not in the commands below) ----
CMS_CONDITIONS="${CMS_CONDITIONS:-auto:phase1_2024_realistic}"   # production GT: 140X_mcRun3_2024_realistic_v26
CMS_ERA="${CMS_ERA:-Run3_2024}"
CMS_HLT_MENU="${CMS_HLT_MENU:-2024v14}"
CMS_BEAMSPOT="${CMS_BEAMSPOT:-DBrealistic}"
CMS_NTHREADS="${CMS_NTHREADS:-4}"
CMS_WRITE_AODSIM="${CMS_WRITE_AODSIM:-0}"   # 1: also write the full AODSIM (step3_reco.root, ~0.4 MB/event); NanoAOD only needs MiniAOD
CMS_PFNANO="${CMS_PFNANO:-1}"       # 1: add ALL PF candidates to NanoAOD (cms/HepMCTools/Relabel/python/pfnano_customise.py)
CMS_TIERS="MINIAODSIM"; [ "${CMS_WRITE_AODSIM}" = "1" ] && CMS_TIERS="AODSIM,MINIAODSIM"
NANO_CUSTOMISE=""; [ "${CMS_PFNANO}" = "1" ] && NANO_CUSTOMISE="--customise HepMCTools/Relabel/pfnano_customise.allPF"

cat > "${OUT}/cms_chain.sh" <<INNER
#!/bin/bash
set -eo pipefail
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=${SCRAM_ARCH}
cd "${WORK}/${CMSSW_RELEASE}/src" && eval \$(scram runtime -sh)
# Work on node-local scratch: ROOT/CMSSW opening files on the EOS FUSE mount from inside the container is
# unreliable (xrootd fallback with a broken URL), and local disk is faster. Results are copied to OUT at the end;
# existing step outputs in OUT are copied in first so that finished steps are skipped and used as inputs.
WORK="\${CMS_SCRATCH:-\${TMPDIR:-/tmp}}/cms_\$(basename "$(dirname "${OUT}")")_\$\$"; mkdir -p "\${WORK}"
trap 'cp -f "\${WORK}"/step*.root "\${WORK}"/*_cfg.py "\${WORK}"/log.* "\${WORK}"/cmd.* "${OUT}/" 2>/dev/null; rm -rf "\${WORK}"' EXIT
cp -f "${OUT}"/step*.root "\${WORK}/" 2>/dev/null || true
cd "\${WORK}"
echo "[cms] \$CMSSW_VERSION (\$SCRAM_ARCH) on \$(hostname), \$(date); scratch \${WORK}"
ln -sf "${HEPMC}" ./events.hepmc          # hepmc_customise.py reads file:events.hepmc from the CWD (plain ASCII read, fine on EOS)
COMMON="--conditions ${CMS_CONDITIONS} --era ${CMS_ERA} --geometry DB:Extended --mc -n ${NEV} --nThreads ${CMS_NTHREADS} --no_exec"

step() { local name=\$1; shift; echo "[cms] === \$name: \$(date) ==="; echo "\$@" > "cmd.\$name"; /usr/bin/time -v "\$@" > "log.\$name" 2>&1 || { echo "[cms] \$name FAILED (see ${OUT}/log.\$name)"; grep -B2 -A12 -m1 "Exception\|Fatal\|Error" "log.\$name" | head -40; exit 1; }; grep -E "Elapsed \(wall|Maximum resident" "log.\$name" | sed "s/^/[cms] \$name /"; }

if [ ! -s step1_gensim.root ]; then
  cmsDriver.py SingleMuPt10_pythia8_cfi -s GEN,SIM \$COMMON --beamspot ${CMS_BEAMSPOT} \\
      --datatier GEN-SIM --eventcontent RAWSIM \\
      --customise HepMCTools/Relabel/hepmc_customise.customise \\
      --customise_commands "process.RandomNumberGeneratorService.g4SimHits.initialSeed = ${SEED}; process.RandomNumberGeneratorService.VtxSmeared.initialSeed = ${SEED}" \\
      --fileout file:step1_gensim.root --python_filename step1_gensim_cfg.py
  step 1_gensim cmsRun step1_gensim_cfg.py
fi
if [ ! -s step2_digiraw.root ]; then
  cmsDriver.py step2 -s DIGI,L1,DIGI2RAW,HLT:${CMS_HLT_MENU} \$COMMON --pileup NoPileUp \\
      --datatier GEN-SIM-RAW --eventcontent RAWSIM \\
      --filein file:step1_gensim.root --fileout file:step2_digiraw.root --python_filename step2_digiraw_cfg.py
  step 2_digiraw cmsRun step2_digiraw_cfg.py
fi
MINI=step3_reco_inMINIAODSIM.root; [ "${CMS_WRITE_AODSIM}" = "1" ] || MINI=step3_reco.root
if [ ! -s \$MINI ]; then
  cmsDriver.py step3 -s RAW2DIGI,L1Reco,RECO,RECOSIM,PAT \$COMMON \\
      --datatier ${CMS_TIERS} --eventcontent ${CMS_TIERS} \\
      --filein file:step2_digiraw.root --fileout file:step3_reco.root --python_filename step3_reco_cfg.py
  step 3_reco cmsRun step3_reco_cfg.py
fi
if [ ! -s step4_nano.root ]; then
  cmsDriver.py step4 -s NANO \$COMMON --scenario pp --datatier NANOAODSIM --eventcontent NANOAODSIM ${NANO_CUSTOMISE} \\
      --filein file:\$MINI --fileout file:step4_nano.root --python_filename step4_nano_cfg.py
  step 4_nano cmsRun step4_nano_cfg.py
fi
echo "[cms] done: \$(date)"; ls -la *.root; cp -f step*.root *_cfg.py log.* cmd.* "${OUT}/"
INNER
chmod +x "${OUT}/cms_chain.sh"
echo "[cms] wrapper=${CMS_CONTAINER_WRAPPER} release=${CMSSW_RELEASE}/${SCRAM_ARCH} input=${HEPMC} out=${OUT}"
exec "${CMS_CONTAINER_WRAPPER}" --no-home --pwd "${OUT}" $(apptainer_binds) -B "${REPO}" -B "${OUT}" -- bash "${OUT}/cms_chain.sh"
