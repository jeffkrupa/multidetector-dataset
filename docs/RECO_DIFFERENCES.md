# ATLAS vs CMS reconstruction differences: measurement programme

Goal: a particle list for each experiment in which an entry **means the same thing** on both sides, with every remaining
difference measured against the common truth rather than assumed. The official collections stay in the files untouched;
the harmonised list is a documented view derived from them (see "Harmonised view" below).

Sample: `Sep19_prod0`, 100k events each of Z→ee, ttbar, QCD dijet (pT-hat > 150 GeV), no pileup, identical generator
events in both detectors (verified per seed by `scripts/compare_events.py`). Fiducial region for everything below:
|η| < 2.5 (tracker acceptance of both experiments). Forward region: not comparable yet (ATLAS sees ~30 % of the truth ΣpT in
2.4 < |η| < 5, CMS ~80 %), parked.

Status legend: **M** measured, **U** understood (mechanism verified), **H** hypothesis only, **T** to do.

## 1. What the two particle lists are

| | ATLAS `JetETMiss{Charged,Neutral}ParticleFlowObjects` | CMS `PFCands` |
|---|---|---|
| entry | a selected track / a topocluster remnant after charged-shower subtraction | one identified, calibrated particle |
| identity | charged / neutral only | h±, e, μ, γ, h0, HF |
| reproduces the experiment's own jets when clustered as stored? | **no**: 50.5 % of jets within 1 % (M) | **yes**: 99.9 % within 1 %; membership table sums to the raw jet pT within 0.25 % (M) |
| what is needed to reproduce the jets | `WeightPFOTool` weights + |z0 sinθ| < 2 mm: 77 % within 1 %, 91 % within 5 %, median 1.000 (M, U) | PUPPI weight (irrelevant without pileup) |

**How ATLAS builds `AnalysisJets` (read from the Athena 25.0.72 configuration, `JetRecConfig/StandardSmallRJets.py`,
`StandardJetConstits.py`, `JetInputConfig.py`, `JetRecTools/JetPFlowSelectionAlg.h`):**

1. `AntiKt4EMPFlow = JetDefinition("AntiKt", 0.4, cst.GPFlow, ...)` — the input is **`GlobalParticleFlowObjects`**, not the
   `JetETMiss` list. The derivation jets (`AntiKt4EMPFlow_deriv`, what PHYSLITE stores) are a clone with the same input.
2. `GlobalParticleFlowObjects` = `JetPFlowSelectionAlg` applied to the `JetETMiss` list. Defaults: remove **charged** objects
   linked to an **LHMedium electron** or a **Medium muon**; **neutral** objects linked to electrons or muons are **kept**
   (`excludeNeutralElectronFE = excludeNeutralMuonFE = false`); photons and taus are not touched by this algorithm.
   (The earlier belief that Global also drops photon- and tau-linked objects was wrong.)
3. Constituent modifiers `["CorrectPFO", "CHS"]`: `CorrectPFOTool` applies the `WeightPFOTool` weight to charged objects
   and the origin correction to neutral ones; `ChargedHadronSubtractionTool` keeps charged objects associated to the
   primary vertex via the track–vertex-association tool (`UseTrackToVertexTool = True` in the default context), not a
   plain |z0 sinθ| cut.
4. anti-kT R = 0.4 → `JetConstitScaleMomentum` (uncalibrated) → calibration → `AnalysisJets.pt`.

`WeightPFOTool` (EM scale): with `EoverP = TracksExpectedEnergyDeposit / E`, w = 1 for pT < 30 GeV, interpolated to `EoverP`
between 30 and 60 GeV, `EoverP` above; minus `EoverP` if `IsInDenseEnvironment`; 0 above 100 GeV. In ttbar 8.4 % of charged
objects are in a dense environment and only 77.8 % of the stored charged pT survives weights + vertex association.

Closure measured (ttbar, jets pT > 30 GeV, |η| < 2, reclustered / stored uncalibrated pT): `JetETMiss` as stored 50 % within
1 %; `JetETMiss` + weights + vertex + origin correction 84 % / 91 % within 1 % / 5 %; **Global + weights 89 % / 97 %**
(100 events, without vertex association or origin correction). Overlap removal removes 0.4 % of objects per event in ttbar;
those few objects are what spoiled the `JetETMiss` closure. CMS `PFCands` reproduce CMS jets to 0.1 %.

**Consequence:** the self-consistent ATLAS particle list (the one ATLAS itself clusters) is `Global{Charged,Neutral}
ParticleFlowObjects`, which our production did not store; the stored `JetETMiss` list differs from it only by the
charged objects linked to Medium electrons and muons (0.4 % of objects in ttbar, more in Z→ee). Test on Z→ee pending.

## 2. Measured so far (fiducial, per event, mean (rms) of reco / truth ΣpT)

| | Z→ee | ttbar | dijet | status |
|---|---|---|---|---|
| CMS total | 0.91 (0.12) | 0.90 (0.05) | 0.90 (0.06) | M |
| ATLAS total, stored list unweighted | 0.80 (0.26) | 0.86 (0.13) | 0.96 (0.13) | M; process dependence H: jet-core over-count + electrons |
| CMS charged / neutral | 1.00 / 0.76 | 1.04 / 0.70 | 1.05 / 0.71 | M |
| ATLAS charged / neutral, unweighted | 0.66 / 1.14 | 0.82 / 0.96 | 0.83 / 1.19 | M |
| charged particles found (ΔR < 0.03), ATLAS / CMS | 76 % / 92 % | 81 % / 93 % | 81 % / 94 % | M; ATLAS TightPrimary, pT > 0.5 GeV (U) |
| charged candidates without truth partner, ATLAS / CMS | 6 % / 15 % | 7 % / 17 % | 7 % / 18 % | M |
| jets pT > 30 GeV with a partner in the other experiment | 99.1 / 99.4 % | 99.3 / 96.6 % | 99.3 / 99.5 % | M; ttbar CMS dip H: lepton-jets |

Electrons (prompt, pT > 10 GeV, 8k events): raw reconstruction CMS 94–95 %, ATLAS 91 %; isolated prompt electrons agree
between processes and detectors within ~1 %; the disagreement is in busy environments and is driven by isolation inside
the CMS cut-based ID (44 % vs 66–70 % for isolation-free IDs) (M). Electron and jet collections are release defaults;
both jet collections contain the leptons as jets (CMS: 99.9 % of medium electrons > 30 GeV are a stored `Jet`) (M).

Known wrong in the comparison scripts (not in the data): ATLAS "MET" is only the soft core term (T: rebuild);
electron selections are not like-for-like (T: matched working points, below).

## 3. Measurement programme

Every item: both detectors, the three processes, binned in truth pT, |η| and local density (truth ΣpT in ΔR < 0.3),
histogrammed per seed and merged (`scripts/compare_particles.py` pattern), numbers exported as JSON next to the plots.

A. **Per-truth-particle response matrix** (the core). For each stable truth particle class — e±, μ±, γ, charged hadron,
   neutral hadron (n, K0L), in-flight decays (K0S, Λ, hyperons) — what it becomes: probability of each reco class (incl.
   "nothing" and "merged into a neighbour"), pT response and resolution, angular resolution. Charged: one-to-one track
   matching (ΔR < 0.03, pT compatibility). Neutral: energy-in-cone bookkeeping, since calorimeters merge.
B. **Unmatched reco objects by class**: secondaries (truth-linked via ATLAS `truthParticleLink` / CMS track info), fakes,
   duplicates, conversions, nuclear-interaction products.
C. **Energy bookkeeping around each truth particle**: reco charged + neutral energy in a cone vs truth, to expose
   double counting (ATLAS dense environment, electrons) and losses (thresholds, CMS neutral −25–30 %).
D. **Identified objects**: electrons, muons, photons, taus — efficiency / misidentification against truth with *matched
   working points* (equal prompt efficiency in bins of pT and truth isolation, derived on Z→ee, frozen, validated on ttbar
   and dijet); ATLAS offers five discrete likelihood points, CMS a continuous score, so CMS is matched to ATLAS.
E. **Jets as a closure test of any particle list**: clustering the list must reproduce the experiment's own jets
   (acceptance test; CMS passes, ATLAS passes at 77 %/1 % with the official weights).
F. **Event-level**: visible ΣpT, charged fraction, MET (after a proper ATLAS rebuild), all with the same definitions.
G. **Forward region** and **pileup**: after the central, no-pileup picture is understood.

## 4. Harmonised view (to be built from A–E, not before)

Common schema per particle: pT, η, φ, mass, charge, class ∈ {charged hadron, electron, muon, photon, neutral hadron /
neutral cluster}, weight, source flags. CMS: `PFCands` × `puppiWeight`. ATLAS candidate recipe: stored flow objects ×
`WeightPFOTool` weight × vertex association, with objects linked (`FE_*Links`) to a *selected* electron / muon / photon
replaced by that object; flagged rather than deleted so the cleaning can be switched off. Acceptance tests: (i) jets
reproduced (E); (ii) Z→ee fiducial response moves from 0.80 (0.26) towards the CMS 0.91 (0.12); (iii) per-class response
matrices (A) of the two detectors as close as the detectors allow, with the residual tabulated here.
Limits already known: ATLAS neutrals carry no γ / hadron identity (only `EM_PROBABILITY`), ATLAS charged selection is
tighter than CMS's; these are reconstruction differences to be quantified, not removed.

## 5. For the next production

Store: jet–constituent links and the jet-input weights (ATLAS), `Global{Charged,Neutral}ParticleFlowObjects`, the
continuous electron likelihood score; keep AOD/MiniAOD for a subset so derivations can be redone; drop generator events
whose final state violates energy–momentum conservation (17 of 300k here; ATLAS truncates their truth record).
