# CMS chain: HepMC2 → GEN-SIM → DIGI-RAW-HLT → AOD/MiniAOD → NanoAOD

Script: `cms/run_cms.sh`. Plugin: `cms/HepMCTools/Relabel/` (built once by `cms/build.sh`).
Everything below was checked on CVMFS / by running on 2026-09-16; "[inferred]" marks the rest.

## Software and environment
| item | value | why |
|---|---|---|
| release | `CMSSW_14_0_25`, `el9_amd64_gcc12` (`el8_amd64_gcc12` + `cmssw-el8` reproduces McM's arch exactly) | last 14_0_X. It is the **only** series still carrying the frozen 2024 HLT menu (`HLT:2024v14`); 15_0_X has `Fake2` for 2024. McM's Run3Summer24 GEN-SIM/DIGI/RECO used 14_0_18–21; MiniAODv6/NanoAODv15 used 15_0_2. One release covers our whole chain. |
| container | `/cvmfs/cms.cern.ch/common/cmssw-el9` wrapper → `/cvmfs/unpacked.cern.ch/registry.hub.docker.com/cmssw/el9:x86_64` | Docker Hub `cmssw/cmssw:*` images stop at 11_1_3 and the open-data images at 10_6_30; CMS distributes software via CVMFS |
| conditions | `auto:phase1_2024_realistic` (resolves to the 2024 MC global tag), era `Run3_2024`, `--beamspot DBrealistic`, `--geometry DB:Extended` | what `runTheMatrix -w upgrade -l 12834.0` (2024 ttbar no-PU relval) uses. Production GTs: `140X_mcRun3_2024_realistic_v26` (GS/DR), `150X_mcRun3_2024_realistic_v2` (MiniAOD/Nano in 15_0). Frontier access from CERN works without any proxy. |

## Feeding HepMC into CMSSW
* **HepMC3 ASCII cannot be read by any release on CVMFS** (2026-09-16): `MCFileSource3` was merged to master on 2026-09-15
  (cms-sw/cmssw PR #51842, no backport). Hence the generator also writes **HepMC2 `IO_GenEvent` ASCII**, read by `MCFileSource`.
* `MCFileSource` publishes its products as `("source","generator")`, but everything downstream expects a module called
  `generator` (`VtxSmeared` wants `generator:unsmeared`, NanoAOD wants `GenEventInfoProduct`, `GenLumiInfoHeader`
  and `GenRunInfoProduct` from it, otherwise PAT/NANO crash in the Rivet/HTXS modules). The 40-line
  `HepMCSourceRelabeler` EDProducer (`cms/HepMCTools/Relabel/plugins`) does exactly that, and
  `hepmc_customise.py` swaps it in for the dummy Pythia fragment cmsDriver insists on (`SingleMuPt10_pythia8_cfi`).
* `MCFileSource` needs `firstLuminosityBlockForEachRun = cms.untracked.VLuminosityBlockID()` (no fillDescriptions).
* Units: no conversion is done; the file must be `U GEV MM` (our generator writes that).
* Event numbers: the source assigns run 1, events 1..N sequentially, which equals the HepMC event number because the file
  is written 1-based without gaps. NanoAOD event order is not input order when multithreaded: join on the `event` branch.
* Vertex smearing (`VtxSmeared`, GT beamspot) is applied on top of the HepMC vertex, as in production. `genWeight` = 1.

## Steps (all in one release, no pile-up)
1. `GEN,SIM` (Geant4) → `step1_gensim.root` (RAWSIM)
2. `DIGI,L1,DIGI2RAW,HLT:2024v14 --pileup NoPileUp` → `step2_digiraw.root` (no `DATAMIX`/`premix_stage2`, which production uses for pre-mixed pile-up)
3. `RAW2DIGI,L1Reco,RECO,RECOSIM,PAT` → `step3_reco.root` (AODSIM) + `step3_reco_inMINIAODSIM.root`
4. `NANO` with `--eventcontent NANOAODSIM` (flat tree; McM's `NANOEDMAODSIM` is EDM and not uproot-flat) → `step4_nano.root`.
   `CMS_PFNANO=1` (default) adds `--customise PhysicsTools/NanoAOD/custom_btv_cff.BTVCustomNanoAOD`, the in-release PFNano,
   giving `PFCands_*` / `JetPFCands_*` on top of the standard `Jet_*`, `FatJet_*`, `Muon_*`, `Electron_*`, `GenPart_*`, `HLT_*` branches.

## Cost (measured, 3 Z→ee events, `--nThreads 4`, lxplus, warm CVMFS cache; `log.<step>` has the full `/usr/bin/time -v`)
| step | wall | max RSS |
|---|---|---|
| GEN-SIM (Geant4) | 58 s | 1.4 GB |
| DIGI+L1+HLT:2024v14 | 41 s | 3.7 GB |
| RECO+PAT (AOD+MiniAOD) | 82 s | 3.1 GB |
| NANO (+PF candidates) | 30 s | 1.7 GB |
Mostly initialisation at this size; the research measurement on 5 ttbar events gave ~13 s/event GEN-SIM, 5 s DIGI/HLT,
11 s RECO+PAT, 3 s NANO. Verified read-back: `event` = [1,2,3], `nElectron`, `Jet_pt`, `PuppiMET_pt`, `nPFCands` present.

## Knobs (env vars read by run_cms.sh)
`CMS_CONDITIONS`, `CMS_ERA`, `CMS_HLT_MENU`, `CMS_BEAMSPOT`, `CMS_NTHREADS`, `CMS_PFNANO`, `CMS_WORK_AREA` (where the scram
project area lives; default `cms/build`).

## References
runTheMatrix workflows 12834.0 (2024 ttbar noPU), 13034.0 (with PU), 15634.0 (2024 HLT-on-digi, 15_0_X);
McM public API `https://cms-pdmv-prod.web.cern.ch/mcm/public/restapi/requests/get_setup/<prepid>` (e.g. `TOP-RunIII2024Summer24GS-00001`);
https://github.com/cms-sw/cmssw/pull/51842 (HepMC3 source); https://github.com/cms-jet/PFNano.
