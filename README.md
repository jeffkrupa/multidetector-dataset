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
| 2a. ATLAS | `atlas/run_atlas.sh <hepmc> <dir> <seed> <N>` | `x86_64-almalinux9` ATLAS apptainer image + `asetup Athena,25.0.72` | `EVNT`, `HITS`, `AOD`, one `DAOD_PHYSLITE` with constituent-level content |
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

`python3 scripts/compare_events.py output/<sample>` (after `source gen/env_lcg.sh`) runs seven checks and exits non-zero if any fails:

| check | what it proves |
|---|---|
| `hepmc-inputs` | the HepMC3 file read by ATLAS and the HepMC2 file read by CMS contain identical particles (order-independent, 1e-6). If the job did not keep the HepMC text, the events are regenerated from the saved card + seed and the sha256 recorded by the job must match (proves the regenerated file is what both experiments read) |
| `atlas-truth` | every truth particle ATLAS kept in PHYSLITE (`TruthElectrons`, `TruthPhotons`, `TruthNeutrinos`, `TruthBosons…`, `TruthTop`…) exists in the HepMC file |
| `cms-truth` | every `GenPart` CMS kept in NanoAOD exists in the HepMC file |
| `atlas-fulltruth` | every HepMC particle exists in the ATLAS `TruthParticles` (full record, from `DAOD_FTAG1`) |
| `cms-gencands` | every HepMC stable particle with \|η\|<5.9 exists in CMS `GenCands` (all packedGenParticles; CMS stores \|y\|<6) |
| `truth-xcheck` | stable truth e/μ stored by both experiments agree 1:1 |
| `event-join` | N generated = N ATLAS = N CMS = N joined on event number |

It then prints a per-event reco table and writes `<sample>/compare/compare.png` (truth identity, ΔR-matched electron and
jet pT ATLAS vs CMS, MET, jet multiplicity, electron response vs HepMC truth), `compare_pflow.png` (particle-flow level:
charged/neutral candidate multiplicity and ΣpT ATLAS vs CMS in a common fiducial region, multiplicity vs HepMC stable
charged particles, charged pT spectra) and `event_<n>.png` η–φ displays.
Caveats: ATLAS MET is the vector sum of the PHYSLITE core terms (final MET needs METMaker in Athena); PHYSLITE jets are not
overlap-removed against electrons, NanoAOD jets are not either, so a Z→ee event shows two "jets" in both.

Verified result (3 Z→ee events, 2026-09-16, lxplus), all checks PASS:

```
event | ATLAS n_ele lead e pt                jets pt (GeV) | CMS n_ele lead e pt                jets pt (GeV)
    1 |           2      47.2                 [51.0, 46.6] |         2      46.1                 [50.7, 46.3]
    2 |           1      13.8                 [21.9, 18.3] |         1      12.9                       [15.3]
    3 |           0      None                 [79.5, 50.1] |         0      None                 [69.3, 52.9]
```

Also verified on ttbar: 2 events (seed 7) and 20 events (seed 11), all checks PASS; 97 ΔR-matched jets lie on the
ATLAS-vs-CMS identity line (`output/ttbar_13p6TeV_seed11_n20/compare/compare.png`).

Particle-flow level (`compare_pflow.png`): the experiments partition the same energy differently. Within |η|<2.5 and
relative to the visible truth ΣpT (200 ttbar events): ATLAS charged 0.51 + neutral 0.42 = 0.93, CMS charged 0.65 + neutral
0.28 = 0.94, per-event totals correlated at 0.96. ATLAS leaves the energy of high-pT and dense-core tracks in the
calorimeter (neutral FlowElements), CMS assigns it to charged hadrons; CMS also keeps ~25% more soft (<2 GeV) charged
candidates. Jets sum both categories, so they agree while candidate counts do not. Per ttbar event (measured, 20 events): ~400 CMS PF candidates; ~65 charged + ~55 neutral ATLAS Global FlowElements (FTAG1, ~175 kB/event incl. tracks, clusters, truth); the ATLAS
charged multiplicity sits ~20% below the HepMC stable-charged count in acceptance (PFlow track selection), CMS ~10% above
(secondaries, split tracks).

Wall time, 20 ttbar events on lxplus, both experiments with 4 threads (the default):

| step | ATLAS (4 thr) | ATLAS (1 thr) | CMS (4 thr) |
|---|---|---|---|
| HepMC in | 1:32 | 0:58 | (in GEN-SIM) |
| Geant4 | 7:19 | 17:26 | 1:39 (GEN-SIM) |
| digi/HLT + reco | 3:04 | 3:05 | 0:48 + 1:04 |
| analysis format(s) | 2:50 | 4:40 | 0:35 |
| **total** | **~15 min** | **~26 min** | **~4 min** |

ATLAS Geant4 CPU is the same in both modes (~1010 s), so the threading is efficient and the reconstructed output is
identical event by event. Reco and derivation are start-up dominated at this size (1–2 min each). Per-event cost at scale is
dominated by Geant4 (~50 s CPU/event ATLAS, ~15 s CPU/event CMS for ttbar).

## Output formats

| | ATLAS | CMS |
|---|---|---|
| **one analysis file, uproot-readable** | `DAOD_PHYSLITE.out.pool.root` (`CollectionTree`): PHYSLITE analysis objects (`Analysis*AuxDyn.*`, MET terms, b-tagging, truth summaries) **plus**, added by `atlas/python/MultiDetPHYSLITE.py`: all `JetETMiss` charged/neutral FlowElements *before* e/γ/μ/τ overlap removal (the CMS-PF analogue), the full `TruthParticles`/`TruthVertices` record, **all** `InDetTrackParticles` (PHYSLITE's track thinning removed) and `CaloCalTopoClusters`. These added containers are read as `<Container>Aux./<Container>Aux.<var>`. ~110 kB/event for ttbar. | `step4_nano.root` (`Events`): NanoAOD with `PFCands_*` for **every** packed PF candidate (pt, eta, phi, mass, pdgId, charge, puppiWeight, track quality, d0/dz), `JetPFCands_*`, `GenCands_*` | `step4_nano.root`: `PFCands_*` for **every** packed PF candidate (pt, eta, phi, mass, pdgId, charge, puppiWeight, track quality, d0/dz), `JetPFCands_*` jet↔candidate index table, `GenCands_*` |
| full reconstruction kept | `AOD.pool.root` | `step3_reco.root` (AODSIM), `step3_reco_inMINIAODSIM.root` |
| detector-level | `HITS.pool.root` | `step1_gensim.root`, `step2_digiraw.root` (RAW) |
| generator-level | `EVNT.pool.root` | `gen/events.hepmc3` / `.hepmc2` (shared) |

## Known small effects (measured on 1000 ttbar events)

* **Quarkonium → gluons events (2 per 1000)**: when Pythia decays an Υ or ψ(2S) to three gluons, Athena's `FixHepMC`
  strips the parton daughters, `TestHepMC` then rejects the event for a decayed particle without a decay vertex, and
  Gen_tf runs off the end of the input file and fails. The generator now drops events containing a decayed particle
  with parton-only (or no) daughters (counted as `dropped_malformed` in `events.json`), so neither experiment sees
  them; every other event is byte-identical to before. The bias is negligible and the same for both experiments.
* **ATLAS PDG-mass adjustment (~1e-5 of particles)**: when storing the truth record ATLAS resets the mass of broad
  resonances generated off-shell (e.g. an a1(1260) at 1.257 GeV → 1.230 GeV) and rescales that particle and its decay
  products by O(0.1–0.2%). The `atlas-fulltruth` check reports how many particles matched only at <0.5% for this reason.

## Limitations / next steps

* No pile-up in either experiment (deliberate; enabling it needs pre-mixed/minimum-bias inputs and is per-experiment work).
* CMS reads HepMC2 because no CMSSW release on CVMFS reads HepMC3 yet (`MCFileSource3` was merged to master on 2026-09-15).
  The two files are written from the same in-memory event, so content is identical.
* The CMS global tag is `auto:phase1_2024_realistic`; swap in the production tag (`140X_mcRun3_2024_realistic_v26`) via `CMS_CONDITIONS` if bit-level agreement with Run3Summer24 matters.
* Scale-out: each `run_all.sh` invocation is one seed; run many seeds as independent HTCondor jobs (the scripts are self-contained and only need CVMFS + apptainer).

## Running at scale (HTCondor at CERN)

```bash
condor_submit condor/submit.sub                                           # 10 jobs x 500 ttbar events, seeds 1000..1009
condor_submit process=zee_13p6TeV nevents=500 njobs=20 seed0=2000 condor/submit.sub
```

Each job (`condor/job.sh`) ships the repo as a tarball, runs gen + ATLAS + CMS on the node's scratch disk (2 cores by
default, `cpus=1` for single-threaded), and copies results to `output/condor/<process>_n<N>/seed<seed>/` on EOS.
What is kept is configurable with `keep=`:

| `keep=` | files copied back | size / event |
|---|---|---|
| `min` (default) | the ATLAS `DAOD_PHYSLITE`, the CMS NanoAOD, generator provenance (card, seed, cross-section, checksums), all logs/commands/configs | ~0.15 MB |
| `std` | + ATLAS `AOD`, `EVNT`; CMS MiniAOD (re-derive / re-dump later) | ~0.6 MB |
| `hepmc=1` (any level) | + the generated events as gzipped HepMC3 text. Off by default: the full generator record is inside the ATLAS file and the generator is deterministic, so any seed's HepMC can be regenerated in seconds with `gen/run_gen.sh` | +65 kB |
| `all` | + ATLAS `HITS`; CMS GEN-SIM, RAW | ~4 MB |

Content switches: `ATLAS_EXTRA_CONTENT="pflow truth tracks clusters"` (default) selects what is added to the single
PHYSLITE (`""` = plain PHYSLITE); `ATLAS_DERIV_FORMATS="PHYSLITE FTAG1"` adds the official flavour-tagging format;
`ATLAS_PFLOW_DUMP=1` adds the flat `pflow.root`; `CMS_WRITE_AODSIM=1` re-enables the full AODSIM output (off by default:
NanoAOD only needs MiniAOD); `CMS_PFNANO=0` drops PF candidates from NanoAOD.
Budget: ~50 s CPU/event ATLAS + ~15 s CPU/event CMS, so 500 ttbar events ≈ 2.5 h on 8 cores; use `+JobFlavour="workday"`
for ≤200 events/job. Then run the checks and plots over all seeds at once:
`python3 scripts/compare_events.py output/condor/<process>_n<N>/seed*` (plots go to `output/condor/<process>_n<N>/compare_all/`).

## Running at scale (Slurm or another scheduler)

`condor/job.sh` is scheduler-agnostic: it takes its scratch directory and core count from HTCondor, Slurm
(`SLURM_TMPDIR`, `SLURM_CPUS_PER_TASK`) or the shell, gets the repository from a tarball (`REPO_TARBALL`) or a visible
checkout (`REPO`), and delivers results to any writable filesystem (or EOS over xrootd). `slurm/submit.sbatch` is a
job-array front end:

```bash
condor/make_tarball.sh
PROCESS=ttbar_13p6TeV NEVENTS=100 SEED0=3000 KEEP=min sbatch --array=0-19 --cpus-per-task=1 slurm/submit.sbatch
```

Submit from the repo root (Slurm runs a spooled copy of the script, so the repo is taken from `SLURM_SUBMIT_DIR`, or
`REPO_DIR`). Nodes need x86_64, CVMFS, user namespaces and outbound access to the ATLAS/CMS Frontier conditions servers;
the host OS does not matter (on a non-EL9 host the generator step also runs inside the EL9 image, and `apptainer` is
taken from `/cvmfs/oasis.opensciencegrid.org` if the site has none). No CERN account is involved; outputs go to
`FINAL_BASE` (default `output/slurm/<process>_n<N>/seed<seed>/`). Each task works in a directory of its own under
`SCRATCH_BASE` (default: the site's `SLURM_TMPDIR`, else `$PSCRATCH`, else `TMPDIR`) and removes it after delivery.

NERSC Perlmutter (SLES 15, no system apptainer, no node-local disk), verified with 10 x 100 Z→ee events, ~1h10 per job
on the 4 hardware threads the `shared` QOS gives a 6 GB request:

```bash
FINAL_BASE=$CFS/<project>/<user>/multidetector/<tag> PROCESS=zee_13p6TeV NEVENTS=100 SEED0=1000 \
  sbatch -A <project> -C cpu -q shared --time=02:00:00 --array=0-9 slurm/submit.sbatch
```

## Requirements

* x86_64 host with `/cvmfs/{sft,atlas,atlas-condb,cms}.cern.ch` and `apptainer` or user namespaces (lxplus works out
  of the box; elsewhere the CVMFS `apptainer` is used).
* Network access to the ATLAS Frontier conditions service and the CMS conditions (Frontier / CVMFS snapshot).
  No grid certificate is needed.
* Disk: keep `output/` on EOS or local scratch; the container steps bind `/eos`, `/tmp` and the Kerberos cache
  (`/run/user/<uid>`) so EOS FUSE stays writable inside the container.
