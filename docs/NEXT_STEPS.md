# Next steps: validation framework, new processes, CMS tracks/clusters

Handoff brief, 2026-10-06. Read together with `docs/RECO_DIFFERENCES.md` (measurement programme A–G, harmonised view,
"for the next production"), which this brief implements rather than restates. State of the repo: `main` is pushed to
`github.com/jeffkrupa/multidetector-dataset`. One untracked draft exists, `scripts/validate_dataset.py` (single-layer
validator); it is to be restructured, not committed as is.

## What exists and works

* One Pythia8 sample (card + seed) → ATLAS (Athena 25.0.72, full Geant4, MC23 no pile-up) and CMS (CMSSW_14_0_25, 2024
  conditions, no pile-up), joined event by event. 3000 ttbar events validated (`output/condor/ttbar_13p6TeV_n100/seed{1000-1009,2000-2019}`).
* Per job the kept files are: ATLAS `DAOD_PHYSLITE.out.pool.root` (analysis objects + `JetETMiss*ParticleFlowObjects`,
  all `InDetTrackParticles` with truth links, `CaloCalTopoClusters`, full `TruthParticles`; added by
  `atlas/python/MultiDetPHYSLITE.py`, read as `<Container>Aux./<Container>Aux.<var>`), CMS `step4_nano.root` (NanoAOD +
  every packed PF candidate `PFCands_*`, `JetPFCands_*`, `GenCands_*`), generator provenance (card, seed, sha256; HepMC
  text itself is not kept and is regenerated on demand, see `compare_events.hepmc_files`).
* `scripts/compare_events.py` (seven PASS/FAIL identity checks, object-level marginals and efficiencies, multi-seed
  merge), `scripts/compare_particles.py` (particle-level marginals of every PF candidate vs truth, per-seed histograms
  merged with `scripts/compare.sh`), `scripts/studies/` (jet closure, electron ID scans behind RECO_DIFFERENCES.md).
* Batch: `condor/job.sh` (scheduler-agnostic), `condor/submit.sub` (EosSubmit schedd, `cpus=1 keep=min` defaults),
  `slurm/submit.sbatch` + `slurm/compare.sbatch` (runs on non-EL9 Slurm sites such as NERSC Perlmutter: generator step
  and comparison run inside the EL9 image, CVMFS apptainer fallback). Measured: ttbar 100 events/job on 1 core ≈ 64 min
  (ATLAS 51, CMS 13) ≈ 0.011 core-h/event.
* ATLAS Gen_tf already keeps every event CMS takes: `atlas/jobconfig/999002` switches off TestHepMC's vertex-displacement
  and energy/momentum-imbalance rejections (the rest stay on).
* Known from RECO_DIFFERENCES.md: ATLAS jets are clustered from `Global{Charged,Neutral}ParticleFlowObjects` with
  `WeightPFOTool` weights and a |z0 sinθ| < 2 mm vertex association; the stored `JetETMiss*` list alone does not
  reproduce ATLAS jets (50 % within 1 %), with the weights it does (77 % / 91 % within 1 % / 5 %).
* Read `README.md`, `docs/ATLAS.md`, `docs/CMS.md` first; pitfalls that cost time are recorded there (EOS `mkdir -p`,
  Kerberos bind for apptainer, `asetup` vs `set -e`, CMSSW on EOS FUSE, quarkonium→ggg events, ATLAS PDG-mass adjustment,
  FTAG1 Global FlowElements being post-overlap-removal).

## 1. Per-event content to add (build first; every new sample should be produced with it)

This merges RECO_DIFFERENCES.md §5 with the tracks/clusters request. Nothing here changes physics; it is what the
candles and the harmonised view (RECO_DIFFERENCES.md §4) need to read back.

ATLAS (`atlas/python/MultiDetPHYSLITE.py`, one `CONTENT` group each, all on by default):
* `Global{Charged,Neutral}ParticleFlowObjects` next to the `JetETMiss*` ones (jets are built from the Global list).
* Jet–constituent links on the PHYSLITE jets (`AnalysisJetsAuxDyn.constituentLinks`, plus `constituentWeights` if
  present) and the jet-input weights per flow object (`WeightPFOTool` output decoration; find its decoration name in the
  Athena jet config, `JetRecConfig`, and add it to the FlowElement variable list; if it is not persisted, run the tool in
  the derivation job and decorate). Verify the links resolve into the stored Global containers.
* The continuous electron likelihood score (`DFCommonElectronsLHLoose`/`Medium`/`Tight` flags are there; add the
  `LHValue` decoration) so CMS working points can be matched to ATLAS (programme item D).
* `InDetTrackParticles` (all, with truth link) and `CaloCalTopoClusters` are already in.

CMS (via the existing plugin package `cms/HepMCTools/Relabel`):
* `Track_*` from `generalTracks`: pt, eta, phi, charge, chi2, ndof, nValidHits, nPixelHits, dxy, dz (leading PV),
  quality bits. Keep `generalTracks` in the MiniAOD output (`--customise` adding `keep recoTracks_generalTracks_*_*` to
  `process.MINIAODSIMoutput.outputCommands` in `cms/run_cms.sh` step 3), then a `SimpleFlatTableProducer<reco::Track>`
  instantiation attached in the NANO step (`pfnano_customise.py`). ~100 tracks/event, ~4 kB/event.
* `PFClusterECAL_*`, `PFClusterHCAL_*`, `PFClusterHF_*` from `particleFlowClusterECAL/HCAL/HF`: energy,
  correctedEnergy, eta, phi, layer, time, size. Same mechanism (keep statements + `SimpleFlatTableProducer<reco::PFCluster>`).
  ~7 kB/event. These are per-subdetector 2D clusters; ATLAS topoclusters are 3D EM+HAD objects — each is its own
  experiment's cluster primitive, compared to truth, not to each other.
* No truth link for CMS tracks (needs SIM-level TrackingParticles); the validator matches tracks to generator particles
  by ΔR and pT identically in both experiments and uses the ATLAS truth link only as a cross-check of that method.
* The continuous electron MVA score is already in NanoAOD (`Electron_mvaIso`, `Electron_mvaNoIso`).

Generator (`gen/gen_hepmc.cc`): also drop events whose final state violates energy–momentum conservation (sum of status-1
four-vectors vs beams beyond a tolerance; 17 of 300k seen; ATLAS truncates their truth record). Count them in
`events.json` like `dropped_malformed`.

Retention: produce a subset of each process with `keep=std` (AOD, EVNT, MiniAOD kept) so derivations can be redone.

Test on `output/smoke_zee_n3` (reusable GEN-SIM/RAW files; `cms/run_cms.sh` and `atlas/run_atlas.sh` skip finished
steps). Update `docs/CMS.md`, `docs/ATLAS.md`, README output table; commit; `condor/make_tarball.sh`.

## 2. Generator: new processes and a particle gun

Cards in `config/pythia/` (13.6 TeV, Monash, `ParticleDecays:limitTau0 = on, tau0Max = 10.` in every card):
* `zee_13p6TeV` (exists), `zmumu_13p6TeV`, `wenu_13p6TeV`, `wmunu_13p6TeV` (`WeakSingleBoson:ffbar2W`, `24:onMode`),
  `hgg_13p6TeV` (`HiggsSM:gg2H`, `25:onIfAny = 22`), `gammajet_13p6TeV_pthat{50to200,200}` (`PromptPhoton:all`),
  `qcd_dijet_13p6TeV_pthat{50to150,150to500,500to1500,1500}` (`HardQCD:all`, `PhaseSpace:pTHatMin/Max`),
  `minbias_13p6TeV` (`SoftQCD:nonDiffractive`). One slice = one card = one seed directory.
* Particle gun in `gen/gen_hepmc.cc`: card keys `Gun:pdg`, `Gun:ptMin`, `Gun:ptMax` (flat in log pT), `Gun:etaMax`,
  `Gun:nPerEvent` (default 1). Write a HepMC3 event with a single production vertex at the origin, beam particles as in
  Pythia (status 4) so readers do not choke, the gun particle status 1, run info/weights as for Pythia events.
  Species: e±, μ±, γ, π± → cards `gun_{e,mu,gamma,pi}_13p6TeV` (the "13p6TeV" is only for naming consistency).
* ATLAS Gen_tf for gun events: `atlas/jobconfig/999002` already disables the displacement and imbalance tests; a
  single-particle event additionally fails the beam-energy test. Add `atlas/jobconfig/999003/mc.HepMC_particle_gun.py`
  (copy of 999002 plus `testSeq.TestHepMC.BeamEnergyTest = False`, `CmeDifference = 1e12`; property names verified in
  Athena 25.0.72; DSID 900020 in `/cvmfs/atlas.cern.ch/repo/sw/Generators/MCJobOptions/900xxx/900020/` is the production
  precedent), selected with `ATLAS_GEN_JOBCONFIG=999003` in `atlas/run_atlas.sh`. CMS `MCFileSource` reads gun events
  as-is; check `VtxSmeared` still applies.
* Record hard-process flags at validation time from the regenerated HepMC (status 21–23 particles, W/Z/H decay modes).

## 3. Validator: three layers

Restructure `scripts/validate_dataset.py` (reuse its metric code) into:
* **universal** (fixed criteria): the seven identity checks, files open, one schema per experiment, no NaN/inf, dropped
  fraction < 1 %, wall-time outliers, no duplicate (seed, event), one (release, conditions, commit) tuple.
* **reference-relative**: tier 3/4 metrics (visible-energy ratio/corr, jet match fractions and pT scale, electron and
  muon scale, multiplicities, ATLAS tracking efficiency by pT bin, charged ratios to truth by threshold, energy
  partition) and the RECO_DIFFERENCES.md §2 table (fiducial reco/truth ΣpT per class, charged-particle finding
  efficiency, unmatched fractions, jet partner fractions — produced by `compare_particles.py`; reuse its histograms
  rather than recomputing) compared with `refs/<process>.yaml` (`metric: {value, tol, tol_kind: abs|rel}` plus `source: {commit,
  athena, cmssw, conditions, seeds, date}`). `--make-ref` writes the file from a trusted run; every run appends
  `{date, commit, metric: value}` to `refs/<process>.history.json` so drift is visible.
* **candles**: `candles/<process>.py`, one function per candle returning `{name, value_atlas, value_cms, value_truth,
  pass, criterion}`; plus the hard-process flag fractions vs the card's expectation for every process.
* Output `validation.json` with per-layer summaries and per-seed detail; exit non-zero on any FAIL. Muon handling has
  to be added to the loaders (`AnalysisMuonsAuxDyn.*` with the PHYSLITE quality selection; NanoAOD `Muon_*` with
  `Muon_looseId`), mirroring the electron treatment (CMS loose cut-based ID, see `load_cms`).

Candles and their criteria (sizing rationale: peak precision p needs N_reco ≈ (σ/p)², width to 10 % needs ~200):
| process | candle | pass |
|---|---|---|
| Z→ee, Z→μμ | m_ℓℓ peak, width, lepton match fraction | peaks within 0.3 % between experiments and 0.5 % of truth; width ratio 0.7–1.4; match > 95 % |
| W→eν, W→μν | m_T Jacobian edge; MET response vs truth ν pT | edge within 1 %; response slopes within 5 % (absolute MET not a criterion: ATLAS MET is the core-term sum) |
| tt̄ | ≥4-jet fraction; W_had and top peaks from truth-matched jets | fraction within 0.03; peaks within 1 % between experiments, 2 % of generator masses |
| QCD dijet (4 slices) | response vs truth in 4 pT × 3 η bins; dijet balance | CMS/ATLAS response within 2 % per bin; balance width ratio 0.7–1.4 |
| γ+jet (2 slices) | photon–jet balance; photon match | within 1 %; match > 95 % |
| H→γγ | m_γγ peak, width | as Z |
| min-bias | tracking efficiency 0.5–1, 1–2 GeV; charged multiplicity vs truth | within 0.02 per bin of reference; ratio within 0.03 of reference |
| guns e, μ, γ, π | response and resolution vs (pT, η), PF object and cluster | within 1 % / 10 % of reference per bin |

## 4. Production (budget agreed: ~140 core-hours, ~22 000 events)

| process | events | job size | notes |
|---|---|---|---|
| Z→ee, Z→μμ | 1000 each | 200 | |
| W→eν, W→μν | 2500 each | 200 | |
| tt̄ | 3000, done | | reuse; re-run only if the CMS track/cluster tables must be in the reference sample |
| QCD dijet | 500 per slice × 4 | 100 for the top two slices (Geant4 cost ~3× ttbar) | |
| γ+jet | 1000 + 500 | 200 | |
| H→γγ | 1000 | 200 | |
| min-bias | 2000 | 400 | cheap events; start-up dominates |
| guns | 6400 per species × 4 | 400 | 8 log-pT × 4 η bins × 200 |

`module load lxbatch/eossubmit; condor_submit process=<card> nevents=<N> njobs=<J> seed0=<S> cpus=1 keep=min condor/submit.sub`,
one cluster per card, distinct `seed0` per card. Watchers get killed on lxplus under memory pressure: poll with
`condor_q`, do not rely on background loops. After each cluster: run the validator, write the reference with
`--make-ref`, commit `refs/`.

## Conventions to keep

* Versions only in `config/versions.env`; physics only in cards; per-experiment knobs are env vars documented in
  `docs/*.md`; every batch job records its commit (`condor/make_tarball.sh` before submitting).
* ATLAS numbers are MeV in the file, CMS GeV; the loaders convert to GeV.
* Event numbers are 1-based per seed; join key is (seed, event number).
* Commit messages end with the attribution lines in the existing history.
