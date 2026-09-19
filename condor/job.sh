#!/bin/bash
# Batch job wrapper (HTCondor or Slurm): one (process, seed, N) through gen + ATLAS + CMS on the node's scratch disk,
# then the results to their final directory (shared filesystem / EOS FUSE, or xrootd if only that is reachable).
#
#   condor/job.sh <process> <seed> <nevents> <final_output_dir_on_eos> [keep=min|std|all] [hepmc=0|1]
#
# keep=min (default): only what analysis needs -- ATLAS DAOD_PHYSLITE, CMS NanoAOD, generator provenance (card, seed,
#                     cross-section, checksums), every log, command and config.          (~0.15 MB/event)
# hepmc=1           : also keep the generated events as gzipped HepMC3 text (~65 kB/event; HepMC2 is derivable from it).
#                     Off by default: the full generator record is already inside the ATLAS file, and the generator is
#                     deterministic (same card + seed + version -> identical file), so it can be regenerated in seconds.
# keep=std          : + ATLAS AOD and EVNT, CMS MiniAOD (lets you re-derive / re-dump later) (~0.9 MB/event)
# keep=all          : + ATLAS HITS, CMS GEN-SIM / RAW / AODSIM                                (~4.5 MB/event)
#
# Self-contained: the repository arrives as repo.tar.gz (transfer_input_files) and is unpacked and built on the node, so
# nothing on the worker depends on /eos or /afs being mounted. Results are pushed to EOS through the FUSE mount if it is
# there, otherwise with xrdcp over xrootd (Kerberos ticket shipped by MY.SendCredential).
set -o pipefail
PROC="${1:?process}"; SEED="${2:?seed}"; NEV="${3:?nevents}"; FINAL="${4:?final output dir}"; KEEP="${5:-min}"; HEPMC="${6:-0}"
case "${KEEP}" in 0) KEEP=std;; 1) KEEP=all;; min|std|all) ;; *) echo "[job] keep must be min|std|all"; exit 2;; esac
# scheduler-provided scratch dir and core count: HTCondor, Slurm, or plain shell (TMPDIR, nproc)
SCRATCH="${_CONDOR_SCRATCH_DIR:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}}"; mkdir -p "${SCRATCH}"; cd "${SCRATCH}"
CPUS="${_CONDOR_JOB_CPUS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"
echo "[job] host=$(hostname) cpus=${CPUS} scratch=${SCRATCH} start=$(date)"; klist 2>&1 | head -2

# --- repository: the tarball shipped by HTCondor (cwd), or one named in REPO_TARBALL (Slurm), or a visible checkout
[ -f repo.tar.gz ] || { [ -n "${REPO_TARBALL:-}" ] && cp "${REPO_TARBALL}" repo.tar.gz; }
if [ -f repo.tar.gz ]; then rm -rf repo && mkdir -p repo && tar xzf repo.tar.gz -C repo && REPO="${SCRATCH}/repo"; echo "[job] repo from tarball: $(head -1 repo/PROVENANCE)"
elif [ -n "${REPO}" ] && [ -f "${REPO}/run_all.sh" ]; then echo "[job] repo from ${REPO}"
else echo "[job] no repo.tar.gz and no usable REPO"; exit 1; fi
export TMPDIR="${SCRATCH}" CMS_WORK_AREA="${SCRATCH}/cmswork"
source "${REPO}/scripts/common.sh"
OUT="${SCRATCH}/${PROC}_seed${SEED}_n${NEV}"

# --- run: with < 8 cores ATLAS then CMS sequentially using all cores, else concurrently with 4+4 threads
rc=0
if [ "${CPUS}" -lt 8 ]; then
  export ATLAS_NTHREADS="${CPUS}" CMS_NTHREADS="${CPUS}"
  STEPS="gen atlas" "${REPO}/run_all.sh" "${PROC}" "${SEED}" "${NEV}" "${OUT}" || rc=$?
  STEPS="cms"       "${REPO}/run_all.sh" "${PROC}" "${SEED}" "${NEV}" "${OUT}" || rc=$?
else
  "${REPO}/run_all.sh" "${PROC}" "${SEED}" "${NEV}" "${OUT}" || rc=$?
fi
echo "[job] run_all rc=${rc} at $(date); collecting results"

# --- collect what to keep
KEEPDIR="${SCRATCH}/keep/$(basename "${FINAL}")"; mkdir -p "${KEEPDIR}/gen" "${KEEPDIR}/atlas" "${KEEPDIR}/cms"
cp -f "${OUT}"/gen/events.json "${OUT}"/gen/pythia.cmnd "${OUT}"/gen/SHA256SUMS "${OUT}"/gen/gen.log "${KEEPDIR}/gen/" 2>/dev/null
[ "${HEPMC}" = 1 ] && gzip -c "${OUT}/gen/events.hepmc3" > "${KEEPDIR}/gen/events.hepmc3.gz"
for exp in atlas cms; do
  # always: analysis-level outputs + logs, commands, configs
  cp -f "${OUT}/${exp}"/DAOD_*.pool.root "${OUT}/${exp}"/pflow.root "${OUT}/${exp}"/step4_nano.root \
        "${OUT}/${exp}"/log.* "${OUT}/${exp}"/cmd.* "${OUT}/${exp}"/*_cfg.py "${OUT}/${exp}"/*_chain.sh "${OUT}/${exp}"/jobReport.json "${KEEPDIR}/${exp}/" 2>/dev/null
  [ "${KEEP}" != min ] && cp -f "${OUT}/${exp}"/AOD.pool.root "${OUT}/${exp}"/EVNT.pool.root "${OUT}/${exp}"/step3_reco*.root "${KEEPDIR}/${exp}/" 2>/dev/null
  [ "${KEEP}" = all ]  && cp -f "${OUT}/${exp}"/HITS.pool.root "${OUT}/${exp}"/step[12]_*.root "${KEEPDIR}/${exp}/" 2>/dev/null
done
cp -f "${OUT}"/atlas.log "${OUT}"/cms.log "${REPO}/PROVENANCE" "${KEEPDIR}/" 2>/dev/null
echo "[job] keeping $(du -sh "${KEEPDIR}" | cut -f1)"

# --- deliver: a mounted filesystem (shared FS, EOS FUSE) if the destination can be created, else xrootd to EOS
if mkdirp "${FINAL}" 2>/dev/null && [ -d "${FINAL}" ]; then   # mkdirp: EOS-FUSE-safe (scripts/common.sh)
  cp -rf "${KEEPDIR}/." "${FINAL}/" && echo "[job] copied to ${FINAL}" || rc=$((rc + 100))
else
  XR="root://eosuser.cern.ch/${FINAL}"
  xrdfs eosuser.cern.ch mkdir -p "${FINAL}" >/dev/null 2>&1
  ( cd "${KEEPDIR}" && find . -type f | while read -r f; do xrdcp -f -s "${f}" "${XR}/${f#./}" || echo "[job] xrdcp failed: ${f}"; done ) && echo "[job] copied via xrootd to ${XR}" || rc=$((rc + 200))
fi

# --- clean up: HTCondor wipes its scratch dir, Slurm / shared scratch does not. Only what this job created, and only
# once the results are delivered (rc < 100); the scratch dir itself only if it is this job's own and now empty.
if [ "${rc}" -lt 100 ]; then
  rm -rf "${OUT}" "${SCRATCH}/keep" "${SCRATCH}/repo" "${SCRATCH}/repo.tar.gz" "${CMS_WORK_AREA}"
  [ -n "${SLURM_TMPDIR:-}" ] && { cd /; rmdir "${SCRATCH}" 2>/dev/null; }
else
  echo "[job] delivery failed: results left in ${KEEPDIR}"
fi
echo "[job] done rc=${rc} at $(date)"
exit ${rc}
