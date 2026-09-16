# multidetector-dataset

Reconstruct **the same generated events** with **two full detector simulations**: ATLAS (Athena) and CMS (CMSSW).

```
                     ┌──────────────────────────┐
  Pythia8 8.317 ───▶ │ events.hepmc3 (+ .hepmc2) │ ──┬──▶ ATLAS: Gen_tf → Sim_tf (Geant4) → Reco_tf → Derivation_tf → DAOD_PHYSLITE
  (LCG_110, one      │ one file, one seed        │   │
   card, one seed)   └──────────────────────────┘   └──▶ CMS:   cmsDriver GEN-SIM (Geant4) → DIGI-RAW-HLT → RECO → MiniAOD → NanoAOD
```

Every event carries the same event number in both outputs, so they can be joined event-by-event.

## Quick start (lxplus / any EL9 host with CVMFS + apptainer)

```bash
./gen/build.sh                                   # once: compile the generator against the pinned LCG view
./run_all.sh zee_13p6TeV 12345 100               # process, seed, N events  ->  output/zee_13p6TeV_seed12345_n100/
STEPS="gen cms" ./run_all.sh ttbar_13p6TeV 1 50  # run a subset of steps
```

Each step is also a standalone script (re-runnable, idempotent per step, all commands logged next to the outputs):

| step | script | runs in | output |
|---|---|---|---|
| 1. generate | `gen/run_gen.sh <process> <seed> <N> <dir>` | LCG_110 view from CVMFS | `events.hepmc3`, `events.hepmc2`, `events.json`, `SHA256SUMS` |
| 2a. ATLAS | `atlas/run_atlas.sh <hepmc> <dir> <seed> <N>` | `x86_64-almalinux9` ATLAS apptainer image + `asetup Athena,25.0.72` | `EVNT`, `HITS`, `AOD`, `DAOD_PHYSLITE.*.pool.root` |
| 2b. CMS | `cms/run_cms.sh <hepmc> <dir> <seed> <N>` | `cmssw-el9` apptainer wrapper + `CMSSW_14_0_25` | `GEN-SIM`, `DIGI-RAW`, `AOD`, `MiniAOD`, `NanoAOD` |

All versions are pinned in one place: [`config/versions.env`](config/versions.env).
Physics processes are Pythia8 cards in [`config/pythia/`](config/pythia/) (add a card = add a process).

## Why it is built this way

* **HepMC is the hand-off.** The only way to get *identical* events into two experiments' frameworks is to
  generate once and pass a file. Both Athena (`TruthIO/HepMCReadFromFile`) and CMSSW (`MCFileSource`) read
  HepMC ASCII. Same file, same seed, byte-identical output (`gen/test_reproducible.sh`).
* **Containers via CVMFS, not Docker Hub.** Neither experiment publishes a self-contained image of a current
  full-sim release (ATLAS docker images stop at Athena 23.0.x). What *is* reproducible and pinned is:
  a fixed OS image (ATLAS AlmaLinux9 image / CMS `cmssw/el9`) + a fixed release from CVMFS + fixed
  conditions tags. This is exactly how both experiments run production, so the workflow is also SOTA by construction.
* **Full Geant4 simulation** in both (`FullG4MT_QS` for ATLAS, standard `SIM` step for CMS), Run-3 (13.6 TeV,
  2024-era geometry/conditions), **no pile-up** by default so the two reconstructions see literally the same
  particles. Pile-up can be enabled per experiment later (needs pre-mixed inputs).
* **Long-lived particles are left to Geant4.** The Pythia cards set `ParticleDecays:limitTau0 = on, tau0Max = 10 mm`,
  which is what both experiments' production configs do; ATLAS's `TestHepMC` rejects events otherwise.

Details and pitfalls per experiment: [`docs/ATLAS.md`](docs/ATLAS.md), [`docs/CMS.md`](docs/CMS.md).

## Comparing the two reconstructions

`python3 scripts/compare_events.py output/<sample>` (after `source gen/env_lcg.sh`) runs five checks and exits non-zero if any fails:

| check | what it proves |
|---|---|
| `hepmc-inputs` | the HepMC3 file read by ATLAS and the HepMC2 file read by CMS contain identical particles (order-independent, 1e-6) |
| `atlas-truth` | every truth particle ATLAS kept in PHYSLITE (`TruthElectrons`, `TruthPhotons`, `TruthNeutrinos`, `TruthBosons…`, `TruthTop`…) exists in the HepMC file |
| `cms-truth` | every `GenPart` CMS kept in NanoAOD exists in the HepMC file |
| `truth-xcheck` | stable truth e/μ stored by both experiments agree 1:1 |
| `event-join` | N generated = N ATLAS = N CMS = N joined on event number |

It then prints a per-event reco table and writes `<sample>/compare/compare.png` (truth identity, ΔR-matched electron and
jet pT ATLAS vs CMS, MET, jet multiplicity, electron response vs HepMC truth) and `event_<n>.png` η–φ displays.
Caveats: ATLAS MET is the vector sum of the PHYSLITE core terms (final MET needs METMaker in Athena); PHYSLITE jets are not
overlap-removed against electrons, NanoAOD jets are not either, so a Z→ee event shows two "jets" in both.

Verified result (3 Z→ee events, 2026-09-16, lxplus), all checks PASS:

```
event | ATLAS n_ele lead e pt                jets pt (GeV) | CMS n_ele lead e pt                jets pt (GeV)
    1 |           2      47.2                 [51.0, 46.6] |         2      46.1                 [50.7, 46.3]
    2 |           1      13.8                 [21.9, 18.3] |         1      12.9                       [15.3]
    3 |           0      None                 [79.5, 50.1] |         0      None                 [69.3, 52.9]
```

Wall time for 3 events including all framework start-up: ATLAS ~8.5 min (1 thread), CMS ~3.5 min (4 threads).
Per-event cost at scale is dominated by Geant4 (tens of seconds per event in each experiment).

## Output formats

| | ATLAS | CMS |
|---|---|---|
| analysis-level, uproot-readable | `DAOD_PHYSLITE.out.pool.root` (`CollectionTree`) | `step4_nano.root` (`Events`), incl. `PFCands_*` |
| full reconstruction kept | `AOD.pool.root` | `step3_reco.root` (AODSIM), `step3_reco_inMINIAODSIM.root` |
| detector-level | `HITS.pool.root` | `step1_gensim.root`, `step2_digiraw.root` (RAW) |
| generator-level | `EVNT.pool.root` | `gen/events.hepmc3` / `.hepmc2` (shared) |

## Limitations / next steps

* No pile-up in either experiment (deliberate; enabling it needs pre-mixed/minimum-bias inputs and is per-experiment work).
* CMS reads HepMC2 because no CMSSW release on CVMFS reads HepMC3 yet (`MCFileSource3` was merged to master on 2026-09-15).
  The two files are written from the same in-memory event, so content is identical.
* The CMS global tag is `auto:phase1_2024_realistic`; swap in the production tag (`140X_mcRun3_2024_realistic_v26`) via `CMS_CONDITIONS` if bit-level agreement with Run3Summer24 matters.
* Scale-out: each `run_all.sh` invocation is one seed; run many seeds as independent HTCondor jobs (the scripts are self-contained and only need CVMFS + apptainer).

## Requirements

* EL9 host with `/cvmfs/{sft,atlas,atlas-condb,cms}.cern.ch` and `apptainer` (lxplus works out of the box).
* Network access to the ATLAS Frontier conditions service and the CMS conditions (Frontier / CVMFS snapshot).
  No grid certificate is needed.
* Disk: keep `output/` on EOS or local scratch; the container steps bind `/eos`, `/tmp` and the Kerberos cache
  (`/run/user/<uid>`) so EOS FUSE stays writable inside the container.
