# ATLAS chain: HepMC → EVNT → HITS → AOD → DAOD_PHYSLITE

Script: `atlas/run_atlas.sh`. Job option: `atlas/jobconfig/999002/mc.HepMC_external_events.py`.
Everything below was checked against Athena 25.0.72 on CVMFS (2026-09-16); "[inferred]" marks what was
composed from release code/tests rather than executed.

## Software and environment
| item | value | why |
|---|---|---|
| release | `asetup Athena,25.0.72` from `/cvmfs/atlas.cern.ch/repo/sw/software/25.0` | newest 25.0 on CVMFS (2026-09-11); one release covers Gen, Sim, Reco, Derivation |
| container | `/cvmfs/atlas.cern.ch/repo/containers/fs/singularity/x86_64-almalinux9` (what `setupATLAS -c el9` uses) | no public self-contained Athena ≥ 24 docker image exists (Docker Hub `atlas/athena` stops at 23.0.27); ATLAS software is distributed via CVMFS |
| conditions | Frontier (`atlasfrontier-ai.cern.ch`, set by asetup) + POOL payloads on `/cvmfs/atlas-condb.cern.ch` | the sqlite DBRelease on CVMFS predates MC23 tags, so Frontier is required. No grid proxy / X509 needed from CERN. Off-site you may need a Frontier squid. |
| geometry / conditions tag | `ATLAS-R3S-2021-03-02-00` / `OFLCOND-MC23-SDR-RUN3-11-02` | release defaults for Run-3 MC (`AthenaConfiguration/TestDefaults.py`) |

## Steps
1. **HepMC → EVNT** (`Gen_tf.py`, legacy job-option path). The reader is `TruthIO/HepMCReadFromFile`
   (`HepMC3::deduce_reader`, so HepMC2 or HepMC3 ASCII both work; units converted to MeV/mm).
   Pitfalls: the input must sit in the run directory and be named `*.events`
   (Gen_tf globs `*<name>.*ev*ts`); job config dir must be a 6-digit DSID with `mc.HepMC_<x>.py`;
   `--ecmEnergy` must match the beam energy in the file (TestHepMC check); only the non-`--CA` Gen_tf supports HepMC input.
2. **EVNT → HITS** (`Sim_tf.py --CA --simulator FullG4MT_QS`, preInclude `Campaigns.MC23eSimulationMultipleIoV`,
   postInclude `PyJobTransforms.UseFrontier`). Pattern from `Simulation/Tests/ISF_Validation/test/test_RUN3_FullG4MT_QS_ttbar.sh`.
3. **HITS → AOD** (`Reco_tf.py --CA`, preInclude `Campaigns.MC23eNoPileUp`) [inferred from `Campaigns` + `RecJobTransformTests/test_mc23a_13p6TeV.sh`].
   No pile-up, no trigger simulation. Production-style pile-up would need presampled `RDO_BKG` inputs and `--steering doOverlay doRDO_TRIG doTRIGtoALL`.
4. **AOD → DAOD_PHYSLITE** (`Derivation_tf.py --CA --formats PHYSLITE`). PHYSLITE is flat enough for uproot:
   `uproot.open("DAOD_PHYSLITE.out.pool.root:CollectionTree")["AnalysisJetsAuxDyn.pt"]`.

## Cost (measured, 3 Z→ee events, single thread, lxplus, warm CVMFS cache; `log.<step>` has the full `/usr/bin/time -v`)
| step | wall | max RSS | note |
|---|---|---|---|
| Gen_tf (HepMC→EVNT) | 59 s | 1.4 GB | almost all Athena start-up |
| Sim_tf FullG4MT_QS  | 2:32 | 2.9 GB | ~1 min start-up, then O(30 s)/event for Z→ee; literature: 90–185 s/event for ttbar |
| Reco_tf (HITS→AOD)  | 2:42 | 4.4 GB | no pile-up, no trigger |
| Derivation_tf PHYSLITE | 2:10 | 4.3 GB | mostly start-up |
Per-event cost at scale is dominated by Sim; start-up (~1–2 min per transform) amortises over 100s of events per job.
Verified read-back: `EventInfoAuxDyn.eventNumber` = [1,2,3] (same as the HepMC event numbers), reconstructed
`AnalysisElectronsAuxDyn.pt` present.

## Knobs (env vars read by run_atlas.sh)
`ATLAS_GEOMETRY`, `ATLAS_CONDITIONS`, `ATLAS_SIMULATOR` (e.g. `ATLFAST3MT_QS` for fast sim), `ATLAS_SIM_PREINCLUDE`,
`ATLAS_RECO_PREINCLUDE`, `ATLAS_ECM_GEV`, `ATLAS_NTHREADS`, `ATLAS_DERIV_FORMATS` (e.g. `PHYS PHYSLITE`).
