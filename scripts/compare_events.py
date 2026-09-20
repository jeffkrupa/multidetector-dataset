#!/usr/bin/env python3
"""Compare the ATLAS and CMS outputs of one sample, event by event, at generator and reconstruction level.

    source gen/env_lcg.sh                      # needs uproot, awkward, numpy, matplotlib (all in the LCG view)
    python3 scripts/compare_events.py output/<sample>                       # one sample
    python3 scripts/compare_events.py output/condor/ttbar_13p6TeV_n100/seed*   # merge many seeds (plots -> .../compare_all/)
  at scale (everything above is held in memory, fine up to ~10k events): one pass per seed, then merge; scripts/compare.sh
  runs this script in the right environment on any host and slurm/compare.sbatch does both steps as batch jobs
    python3 scripts/compare_events.py --no-plots --summary sum/seed1000.pkl <dir>/seed1000      # checks + slim summary
    python3 scripts/compare_events.py --merge --out <plotdir> sum/*.pkl                         # check tally + plots

Checks (each prints PASS/FAIL and the script exits non-zero on any FAIL):
  1. hepmc-inputs   the .hepmc2 file given to CMS and the .hepmc3 file given to ATLAS contain identical particles
  2. atlas-truth    every truth particle ATLAS stored (PHYSLITE Truth* containers) exists in the HepMC file
  3. cms-truth      every gen particle CMS stored (NanoAOD GenPart) exists in the HepMC file
  2b. atlas-fulltruth  every HepMC particle exists in the ATLAS TruthParticles (full record, from DAOD_FTAG1 or pflow.root)
  3b. cms-gencands     every HepMC stable particle with |eta|<5.9 exists in the CMS NanoAOD GenCands (all packedGenParticles, |y|<6)
  4. truth-xcheck   stable truth leptons/photons stored by both experiments agree 1:1 in (pdgId, pt, eta, phi)
  5. event-join     both experiments wrote every generated event, joined on the event number

Then a per-event reco table and plots in <sample>/compare/: compare.png (truth identity, reco vs reco, reco vs truth),
compare_pflow.png (particle-flow candidates: ATLAS Global*ParticleFlowObjects from DAOD_FTAG1 vs CMS PFCands from NanoAOD),
compare_marginals.png (every object of each experiment, matched or not: marginal distributions, reconstruction efficiency
against the common truth, and how often an object of one experiment has a partner in the other) and event_<n>.png displays.
Event numbers: both experiments store the 1-based HepMC event number. NanoAOD may reorder events; nothing assumes order.
"""
import glob, json, math, os, sys
import numpy as np
import awkward as ak
import uproot

# --------------------------------------------------------------------------------------- HepMC (no deps)
def read_hepmc(path):
    """Minimal HepMC2/HepMC3 ASCII reader -> {event_number: [(pdg, status, px, py, pz, e), ...]} in GeV."""
    events, cur, evno, v3 = {}, None, None, None
    import gzip
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        for line in f:
            tag = line[:2]
            if tag == "E ":
                if cur is not None: events[evno] = cur
                cur = []
                v3 = line.split()[1].isdigit() and len(line.split()) == 4  # HepMC3: "E n nvtx npart"
                evno = int(line.split()[1])
            elif tag == "P " and cur is not None:
                s = line.split()
                if v3:   # HepMC3: P id mother pdg px py pz e m status
                    cur.append((int(s[3]), int(s[9]), float(s[4]), float(s[5]), float(s[6]), float(s[7])))
                else:    # HepMC2: P barcode pdg px py pz e m status ...
                    cur.append((int(s[2]), int(s[8]), float(s[3]), float(s[4]), float(s[5]), float(s[6])))
    if cur is not None: events[evno] = cur
    return events

def kin(px, py, pz, e):
    pt = math.hypot(px, py); p = math.sqrt(px * px + py * py + pz * pz)
    eta = 0.5 * math.log((p + pz) / (p - pz)) if p > abs(pz) else (99. if pz > 0 else -99.)
    return pt, eta, math.atan2(py, px)

def match_to_hepmc(cands, hep, rel_tol):
    """cands: list of (pdg, pt, eta, phi). Return list of unmatched candidates (pt above 0.5 GeV only)."""
    index = {}
    for pdg, st, px, py, pz, e in hep:
        index.setdefault(pdg, []).append(kin(px, py, pz, e))
    bad = []
    for pdg, pt, eta, phi in cands:
        if pt < 0.5: continue
        ok = any(abs(p - pt) <= rel_tol * max(pt, 1.) and abs(h - eta) < 0.02 and abs((f - phi + math.pi) % (2 * math.pi) - math.pi) < 0.02
                 for p, h, f in index.get(pdg, []))
        if not ok: bad.append((pdg, round(pt, 2), round(eta, 3), round(phi, 3)))
    return bad

def match_lists(cands, ref, rel_tol):
    """Every (pdg, pt, eta, phi) in cands must exist in ref within tolerance; returns the unmatched ones (pt > 0.5 GeV only)."""
    index = {}
    for pdg, pt, eta, phi in ref: index.setdefault(pdg, []).append((pt, eta, phi))
    bad = []
    for pdg, pt, eta, phi in cands:
        if pt < 0.5: continue
        if not any(abs(p - pt) <= rel_tol * max(pt, 1.) and abs(h - eta) < 0.02 and abs((f - phi + math.pi) % (2 * math.pi) - math.pi) < 0.02
                   for p, h, f in index.get(pdg, [])):
            bad.append((pdg, round(pt, 2), round(eta, 3), round(phi, 3)))
    return bad

# --------------------------------------------------------------------------------------- ATLAS
def load_atlas(d):
    f = glob.glob(os.path.join(d, "atlas", "DAOD_PHYSLITE.*.pool.root"))
    if not f: return None
    t = uproot.open(f[0])["CollectionTree"]; keys = set(t.keys())
    def arr(b): return t[b].array() if b in keys else None
    evn = arr("EventInfoAuxDyn.eventNumber")
    out = {}
    truth_containers = [c for c in ["TruthElectrons", "TruthMuons", "TruthPhotons", "TruthNeutrinos", "TruthTop",
                                    "TruthBosonsWithDecayParticles", "TruthTaus", "TruthBottom"] if f"{c}AuxDyn.pdgId" in keys]
    tr = {c: {v: arr(f"{c}AuxDyn.{v}") for v in ["pdgId", "status", "px", "py", "pz", "e"]} for c in truth_containers}
    reco = {v: arr(v) for v in ["AnalysisElectronsAuxDyn.pt", "AnalysisElectronsAuxDyn.eta", "AnalysisElectronsAuxDyn.phi",
                                "AnalysisMuonsAuxDyn.pt", "AnalysisJetsAuxDyn.pt", "AnalysisJetsAuxDyn.eta", "AnalysisJetsAuxDyn.phi",
                                "MET_Core_AnalysisMETAuxDyn.mpx", "MET_Core_AnalysisMETAuxDyn.mpy", "MET_Core_AnalysisMETAuxDyn.source"]}
    for i in range(len(evn)):
        ev = int(evn[i]); g = lambda v, default=[]: (reco[v][i] if reco[v] is not None else default)
        truth = []   # (pdg, status, pt, eta, phi) in GeV
        for c in truth_containers:
            for j in range(len(tr[c]["pdgId"][i])):
                pt, eta, phi = kin(*(float(tr[c][v][i][j]) / 1000. for v in ["px", "py", "pz", "e"]))
                truth.append((int(tr[c]["pdgId"][i][j]), int(tr[c]["status"][i][j]), pt, eta, phi))
        # PHYSLITE stores MET *terms* (MET_Core_AnalysisMET); the final MET needs METMaker inside Athena.
        # As a framework-free proxy we sum the stored terms vectorially.
        mpx, mpy = g("MET_Core_AnalysisMETAuxDyn.mpx"), g("MET_Core_AnalysisMETAuxDyn.mpy")
        met = math.hypot(float(ak.sum(mpx)), float(ak.sum(mpy))) / 1000. if len(mpx) else None
        out[ev] = dict(truth=truth,
                       ele=[(float(p) / 1000., float(e), float(f)) for p, e, f in zip(g("AnalysisElectronsAuxDyn.pt"), g("AnalysisElectronsAuxDyn.eta"), g("AnalysisElectronsAuxDyn.phi"))],
                       n_mu=len(g("AnalysisMuonsAuxDyn.pt")),
                       jets=[(float(p) / 1000., float(e), float(f)) for p, e, f in zip(g("AnalysisJetsAuxDyn.pt"), g("AnalysisJetsAuxDyn.eta"), g("AnalysisJetsAuxDyn.phi"))],
                       met=met)
    # particle-flow objects + full truth. Preferred source: pflow.root (JetETMiss* FlowElements = all charged/neutral
    # PFlow objects BEFORE e/gamma/mu/tau overlap removal, the analogue of CMS packed PF candidates). Fallback: DAOD_FTAG1
    # Global* FlowElements, which are AFTER overlap removal (tracks matched to any electron/muon/tau *candidate* are gone,
    # so high-pT jet-core tracks are under-counted by ~30%). Truth comes from either.
    ft = glob.glob(os.path.join(d, "atlas", "DAOD_FTAG1.*.pool.root")); pf = os.path.join(d, "atlas", "pflow.root")
    sa = lambda c, v: f"{c}Aux./{c}Aux.{v}"          # split static-aux branch naming used for the added containers
    if sa("JetETMissChargedParticleFlowObjects", "pt") in keys and sa("TruthParticles", "pdgId") in keys:
        a = t.arrays([sa("JetETMissChargedParticleFlowObjects", "pt"), sa("JetETMissChargedParticleFlowObjects", "eta"),
                      sa("JetETMissNeutralParticleFlowObjects", "pt"), sa("JetETMissNeutralParticleFlowObjects", "eta"),
                      sa("TruthParticles", "pdgId"), sa("TruthParticles", "px"), sa("TruthParticles", "py"), sa("TruthParticles", "pz"), sa("TruthParticles", "e")])
        for i in range(len(evn)):
            ev = int(evn[i])
            out[ev]["pf"] = dict(ch=[(float(p) / 1000., float(e)) for p, e in zip(a[sa("JetETMissChargedParticleFlowObjects", "pt")][i], a[sa("JetETMissChargedParticleFlowObjects", "eta")][i])],
                                 ne=[(float(p) / 1000., float(e)) for p, e in zip(a[sa("JetETMissNeutralParticleFlowObjects", "pt")][i], a[sa("JetETMissNeutralParticleFlowObjects", "eta")][i])])
            out[ev]["fulltruth"] = [(int(g), *kin(float(px) / 1000., float(py) / 1000., float(pz) / 1000., float(e) / 1000.)) for g, px, py, pz, e in
                                    zip(a[sa("TruthParticles", "pdgId")][i], a[sa("TruthParticles", "px")][i], a[sa("TruthParticles", "py")][i], a[sa("TruthParticles", "pz")][i], a[sa("TruthParticles", "e")][i])]
            out[ev]["pf_source"] = "PHYSLITE+MultiDet (JetETMiss FE)"
    elif os.path.exists(pf):
        a = uproot.open(pf)["pflow"].arrays(["event", "fe_ch_pt", "fe_ch_eta", "fe_ne_pt", "fe_ne_eta", "tp_pdg", "tp_pt", "tp_eta", "tp_phi"])
        for i in range(len(a)):
            ev = int(a["event"][i])
            if ev in out:
                out[ev]["pf"] = dict(ch=[(float(p), float(e)) for p, e in zip(a["fe_ch_pt"][i], a["fe_ch_eta"][i])],
                                     ne=[(float(p), float(e)) for p, e in zip(a["fe_ne_pt"][i], a["fe_ne_eta"][i])])
                out[ev]["fulltruth"] = [(int(g), float(pt), float(eta), float(phi)) for g, pt, eta, phi in zip(a["tp_pdg"][i], a["tp_pt"][i], a["tp_eta"][i], a["tp_phi"][i])]
                out[ev]["pf_source"] = "pflow.root (JetETMiss FE)"
    elif ft:
        f = uproot.open(ft[0])["CollectionTree"]
        a = f.arrays(["EventInfoAuxDyn.eventNumber", "GlobalChargedParticleFlowObjectsAuxDyn.pt", "GlobalChargedParticleFlowObjectsAuxDyn.eta",
                      "GlobalNeutralParticleFlowObjectsAuxDyn.pt", "GlobalNeutralParticleFlowObjectsAuxDyn.eta",
                      "TruthParticlesAuxDyn.pdgId", "TruthParticlesAuxDyn.px", "TruthParticlesAuxDyn.py", "TruthParticlesAuxDyn.pz", "TruthParticlesAuxDyn.e"])
        for i in range(len(a)):
            ev = int(a["EventInfoAuxDyn.eventNumber"][i])
            if ev not in out: continue
            out[ev]["pf"] = dict(ch=[(float(p) / 1000., float(e)) for p, e in zip(a["GlobalChargedParticleFlowObjectsAuxDyn.pt"][i], a["GlobalChargedParticleFlowObjectsAuxDyn.eta"][i])],
                                 ne=[(float(p) / 1000., float(e)) for p, e in zip(a["GlobalNeutralParticleFlowObjectsAuxDyn.pt"][i], a["GlobalNeutralParticleFlowObjectsAuxDyn.eta"][i])])
            out[ev]["fulltruth"] = [(int(g), *kin(float(px) / 1000., float(py) / 1000., float(pz) / 1000., float(e) / 1000.)) for g, px, py, pz, e in
                                    zip(a["TruthParticlesAuxDyn.pdgId"][i], a["TruthParticlesAuxDyn.px"][i], a["TruthParticlesAuxDyn.py"][i], a["TruthParticlesAuxDyn.pz"][i], a["TruthParticlesAuxDyn.e"][i])]
            out[ev]["pf_source"] = "FTAG1 (Global FE, post overlap removal)"
    return out, truth_containers

# --------------------------------------------------------------------------------------- CMS
def load_cms(d):
    f = os.path.join(d, "cms", "step4_nano.root")
    if not os.path.exists(f): return None
    t = uproot.open(f)["Events"]; keys = set(t.keys())
    want = ["event", "GenPart_pdgId", "GenPart_status", "GenPart_pt", "GenPart_eta", "GenPart_phi",
            "Electron_pt", "Electron_eta", "Electron_phi", "Electron_cutBased", "nMuon", "Jet_pt", "Jet_eta", "Jet_phi", "PuppiMET_pt", "GenMET_pt",
            "PFCands_pt", "PFCands_eta", "PFCands_charge", "GenCands_pdgId", "GenCands_pt", "GenCands_eta", "GenCands_phi"]
    a = t.arrays([w for w in want if w in keys])
    out = {}
    for i in range(len(a)):
        ev = int(a["event"][i])
        truth = [(int(p), int(s), float(pt), float(eta), float(phi)) for p, s, pt, eta, phi in
                 zip(a["GenPart_pdgId"][i], a["GenPart_status"][i], a["GenPart_pt"][i], a["GenPart_eta"][i], a["GenPart_phi"][i])]
        out[ev] = dict(truth=truth,
                       # NanoAOD keeps every electron candidate (no ID); PHYSLITE AnalysisElectrons carry a loose LH ID.
                       # Require the NanoAOD loose cut-based ID (cutBased >= 2) so the two collections are comparable.
                       ele=[(float(p), float(e), float(f)) for p, e, f, idv in zip(a["Electron_pt"][i], a["Electron_eta"][i], a["Electron_phi"][i], a["Electron_cutBased"][i]) if idv >= 2],
                       n_mu=int(a["nMuon"][i]),
                       jets=[(float(p), float(e), float(f)) for p, e, f in zip(a["Jet_pt"][i], a["Jet_eta"][i], a["Jet_phi"][i])],
                       met=float(a["PuppiMET_pt"][i]), genmet=float(a["GenMET_pt"][i]) if "GenMET_pt" in keys else None)
        if "GenCands_pt" in keys:
            out[ev]["gencands"] = [(int(g), float(pt), float(eta), float(phi)) for g, pt, eta, phi in zip(a["GenCands_pdgId"][i], a["GenCands_pt"][i], a["GenCands_eta"][i], a["GenCands_phi"][i])]
        if "PFCands_pt" in keys:
            out[ev]["pf"] = dict(ch=[(float(p), float(e)) for p, e, q in zip(a["PFCands_pt"][i], a["PFCands_eta"][i], a["PFCands_charge"][i]) if q != 0],
                                 ne=[(float(p), float(e)) for p, e, q in zip(a["PFCands_pt"][i], a["PFCands_eta"][i], a["PFCands_charge"][i]) if q == 0])
    return out

# --------------------------------------------------------------------------------------- helpers
def dr(a, b):
    dphi = (a[2] - b[2] + math.pi) % (2 * math.pi) - math.pi
    return math.hypot(a[1] - b[1], dphi)

def match_objects(A, B, maxdr):
    """Greedy ΔR matching of two lists of (pt, eta, phi). Returns list of (a, b)."""
    pairs, used = [], set()
    for x in sorted(A, key=lambda o: -o[0]):
        best = min(((dr(x, y), j) for j, y in enumerate(B) if j not in used), default=None)
        if best and best[0] < maxdr: pairs.append((x, B[best[1]])); used.add(best[1])
    return pairs

RESULTS = []   # (check, passed, detail) of this run, kept for --summary
TITLE = None   # plot title prefix; default: name of the directory above the plot directory
def report(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    RESULTS.append((name, bool(ok), detail)); return ok

# --------------------------------------------------------------------------------------- main
def hepmc_files(d, gen):
    """Return (hepmc3 path, hepmc2 path, note). Batch jobs keep no HepMC text by default (or only gzipped HepMC3), so
    regenerate from the saved card + seed + N with the same generator build (deterministic) and verify the sha256 that
    the job recorded in gen/SHA256SUMS: that proves the regenerated file is byte-identical to what both experiments read."""
    import hashlib, subprocess, tempfile
    g = os.path.join(d, "gen")
    h3 = next((p for p in (g + "/events.hepmc3", g + "/events.hepmc3.gz") if os.path.exists(p)), None)
    h2 = g + "/events.hepmc2" if os.path.exists(g + "/events.hepmc2") else None
    if h3 and h2: return h3, h2, ""
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tmp = tempfile.mkdtemp(prefix="hepmc_regen_", dir=os.environ.get("TMPDIR"))
    if not os.path.exists(repo + "/gen/build/gen_hepmc"): subprocess.run(["bash", repo + "/gen/build.sh"], check=True, capture_output=True)
    subprocess.run([repo + "/gen/build/gen_hepmc", g + "/pythia.cmnd", str(gen["seed"]), str(gen["nevents"]), tmp + "/events"], check=True, capture_output=True)
    recorded = {l.split()[1].split("/")[-1]: l.split()[0] for l in open(g + "/SHA256SUMS")} if os.path.exists(g + "/SHA256SUMS") else {}
    ok = all(hashlib.sha256(open(f"{tmp}/events.{x}", "rb").read()).hexdigest() == recorded.get(f"events.{x}") for x in ("hepmc3", "hepmc2") if f"events.{x}" in recorded)
    note = f"HepMC text not kept; regenerated {gen['nevents']} events from card+seed {gen['seed']}: sha256 {'MATCHES' if ok else 'DOES NOT MATCH'} the job's SHA256SUMS"
    if not ok: raise RuntimeError(note)
    return tmp + "/events.hepmc3", tmp + "/events.hepmc2", note

def load_sample(d):
    """Load one sample directory; keys are (sample_tag, event_number) so several seeds can be merged."""
    d = d.rstrip("/"); tag = os.path.basename(d)
    gen = json.load(open(os.path.join(d, "gen", "events.json")))
    rekey = lambda m: {(tag, k): v for k, v in m.items()} if m is not None else None
    h3p, h2p, note = hepmc_files(d, gen)
    h3 = rekey(read_hepmc(h3p)); h2 = rekey(read_hepmc(h2p))
    if note: print(f"[{tag}] {note}")
    atlas = load_atlas(d); cms = rekey(load_cms(d))
    atlas, containers = (rekey(atlas[0]), atlas[1]) if atlas else (None, [])
    return gen, h3, h2, atlas, cms, containers

def main(dirs, outdir=None, summary=None, plots=True):
    dirs = [d.rstrip("/") for d in dirs]; ok = True
    gen = dict(nevents=0, seeds=[]); h3, h2, atlas, cms = {}, {}, {}, {}; atlas_containers = []; missing = []
    for d in dirs:
        g, a3, a2, at, cm, cont = load_sample(d)
        gen["nevents"] += g["nevents"]; gen["seeds"].append(g["seed"]); gen.setdefault("pythia_version", g["pythia_version"]); gen.setdefault("card", g["card"])
        h3.update(a3); h2.update(a2)
        if at is None or cm is None: missing.append((d, "ATLAS" if at is None else "", "CMS" if cm is None else "")); continue
        atlas.update(at); cms.update(cm); atlas_containers = cont
    print(f"samples: {len(dirs)} dir(s), seeds {gen['seeds'] if len(dirs) <= 12 else str(gen['seeds'][:12]) + '...'}\n"
          f"generator: Pythia {gen['pythia_version']} card {os.path.basename(gen['card'])}, {gen['nevents']} events generated\n")
    for m in missing: print(f"MISSING output in {m[0]}: {m[1]} {m[2]}")

    # 1. HepMC inputs identical
    key = lambda p: (p[0], p[1], p[5], p[4], p[3], p[2])   # order-independent (HepMC2 lists particles per vertex); full precision so near-identical soft partons pair correctly
    same = set(h3) == set(h2) and all(len(h3[e]) == len(h2[e]) and all(
        p[0] == q[0] and p[1] == q[1] and all(abs(p[k] - q[k]) <= 1e-6 * max(1., abs(p[k])) for k in range(2, 6))
        for p, q in zip(sorted(h3[e], key=key), sorted(h2[e], key=key))) for e in h3)
    ok &= report("hepmc-inputs", same, f"{len(h3)} events, {sum(len(v) for v in h3.values())} particles, "
                 f"{sum(sum(1 for p in v if p[1] == 1) for v in h3.values())} stable; HepMC3 (ATLAS) == HepMC2 (CMS)")

    if not atlas or not cms:
        print("no sample with both ATLAS and CMS output"); return 1

    # 2./3. stored truth ⊂ HepMC
    for name, coll, tol in [("atlas-truth", atlas, 5e-3), ("cms-truth", cms, 0.01)]:   # NanoAOD stores pt/eta/phi at ~1e-3 precision; ATLAS: see PDG-mass note below
        nchk, bad = 0, []
        for ev, rec in coll.items():
            c = [(p, pt, eta, phi) for p, s, pt, eta, phi in rec["truth"]]
            nchk += len(c); bad += [(ev, b) for b in match_to_hepmc(c, h3.get(ev, []), tol)]
        src = f"PHYSLITE {'+'.join(atlas_containers)}" if name == "atlas-truth" else "NanoAOD GenPart"
        ok &= report(name, not bad, f"{nchk} stored truth particles ({src}) all found in HepMC" if not bad else f"{len(bad)} unmatched, e.g. {bad[:3]}")

    # 2b/3b. the other direction, full record: every HepMC particle is in ATLAS TruthParticles (pflow.root), and every
    #        HepMC *stable* particle is in CMS GenCands (all packedGenParticles). ATLAS may hold extra Geant4-added
    #        secondaries (decays in flight, conversions), CMS may not, hence subset checks rather than equality.
    if all("fulltruth" in atlas[ev] for ev in atlas):
        # ATLAS resets the mass of broad resonances generated off-shell (e.g. a1(1260)) to the PDG value and rescales
        # that particle and its decay products by O(0.1%) when it stores the truth record. Match at 1e-3 first; what is
        # left must still match at 5e-3 (that residual is the PDG-mass adjustment and is reported, not failed).
        nchk, strict, loose = 0, [], []
        for ev, rec in atlas.items():
            c = [(p, *kin(px, py, pz, e)) for p, st, px, py, pz, e in h3.get(ev, [])]
            nchk += len(c); s1 = match_lists(c, rec["fulltruth"], 1e-3); strict += [(ev, b) for b in s1]
            loose += [(ev, b) for b in match_lists([(b[0], b[1], b[2], b[3]) for b in s1], rec["fulltruth"], 5e-3)]
        # Athena's FixHepMC removes a few intermediate partons per ~10^6 particles "while cleaning decay chains" (seen in events
        # whose Pythia record has an energy-momentum imbalance). Quarks, gluons and diquarks never reach the detector, so
        # that is reported, not failed; anything else missing fails.
        is_parton = lambda g: abs(g) in (1, 2, 3, 4, 5, 6, 21) or (1000 < abs(g) < 10000 and (abs(g) // 10) % 10 == 0)
        partons = [x for x in loose if is_parton(x[1][0])]; loose = [x for x in loose if not is_parton(x[1][0])]
        src = next(iter(atlas.values()))["pf_source"]
        ok &= report("atlas-fulltruth", not loose, f"all {nchk} HepMC particles found in ATLAS TruthParticles ({src}); "
                     f"{len(strict)} matched only at <0.5% (ATLAS PDG-mass adjustment of off-shell resonances)"
                     + (f"; except {len(partons)} intermediate partons removed by ATLAS FixHepMC in events {sorted({x[0][1] for x in partons})}" if partons else "") if not loose
                     else f"{len(loose)} HepMC particles missing from ATLAS truth even at 0.5%, e.g. {loose[:3]}")
    if all("gencands" in cms[ev] for ev in cms):
        nchk, bad = 0, []
        for ev, rec in cms.items():
            # CMS packedGenParticles are produced with maxRapidity = 6 (PATPackedGenParticleProducer); stay inside it
            c = [(p, *kin(px, py, pz, e)) for p, st, px, py, pz, e in h3.get(ev, []) if st == 1 and abs(kin(px, py, pz, e)[1]) < 5.9]
            nchk += len(c); bad += [(ev, b) for b in match_lists(c, rec["gencands"], 0.01)]
        ok &= report("cms-gencands", not bad, f"all {nchk} HepMC stable particles with |eta|<5.9 found in CMS GenCands (packedGenParticles)" if not bad else f"{len(bad)} HepMC stable particles missing from CMS GenCands, e.g. {bad[:3]}")

    # 4. stable leptons/photons stored by both agree 1:1
    n_pairs, n_miss = 0, 0
    for ev in sorted(set(atlas) & set(cms)):
        A = [(pt, eta, phi, p) for p, s, pt, eta, phi in atlas[ev]["truth"] if s == 1 and abs(p) in (11, 13) and pt > 5]
        C = {(p, s) for p, s, pt, eta, phi in cms[ev]["truth"]}
        for pt, eta, phi, p in A:
            n_pairs += 1
            if not any(q == p and s == 1 and abs(pt2 - pt) < 0.01 * pt and dr((pt, eta, phi), (pt2, eta2, phi2)) < 0.01
                       for q, s, pt2, eta2, phi2 in cms[ev]["truth"]): n_miss += 1
    ok &= report("truth-xcheck", n_miss == 0, f"{n_pairs} stable truth e/mu (pT>5 GeV) in ATLAS PHYSLITE, {n_pairs - n_miss} found identically in CMS GenPart")

    # 5. event join
    common = sorted(set(atlas) & set(cms))
    ok &= report("event-join", len(common) == gen["nevents"] == len(atlas) == len(cms) and not missing,
                 f"generated {gen['nevents']}, ATLAS {len(atlas)}, CMS {len(cms)}, joined {len(common)}")

    # per-event reco table
    print(f"\n{'event':>5} | {'ATLAS n_e':>9} {'e1 pT':>6} {'n_jet':>5} {'jet pT (GeV)':>24} {'MET':>6} | {'CMS n_e':>7} {'e1 pT':>6} {'n_jet':>5} {'jet pT (GeV)':>24} {'MET':>6}")
    fmt = lambda js: str([round(j[0], 1) for j in js[:4]])
    for ev in (common if len(common) <= 40 else common[:20]):
        a, c = atlas[ev], cms[ev]
        print(f"{ev[1]:>5} | {len(a['ele']):>9} {(round(a['ele'][0][0],1) if a['ele'] else '-'):>6} {len(a['jets']):>5} {fmt(a['jets']):>24} {(round(a['met'],1) if a['met'] is not None else '-'):>6} | "
              f"{len(c['ele']):>7} {(round(c['ele'][0][0],1) if c['ele'] else '-'):>6} {len(c['jets']):>5} {fmt(c['jets']):>24} {round(c['met'],1):>6}")

    if len(common) > 40: print(f"  ... ({len(common)} events; table truncated)")

    if summary: write_summary(summary, gen, common, atlas, cms, h3)
    if plots:
        outdir = outdir or (os.path.join(dirs[0], "compare") if len(dirs) == 1 else os.path.join(os.path.dirname(dirs[0]), "compare_all")); os.makedirs(outdir, exist_ok=True)
        make_plots(common, atlas, cms, h3, outdir)
        print(f"\nplots: {outdir}/")
    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED"); return 0 if ok else 2

def write_summary(path, gen, common, atlas, cms, hep):
    """Pickle what the plots need and nothing else: the reco objects and stored truth of both experiments (no full truth
    records, particle flow inside the tracker only) and the HepMC stable particles in the tracker acceptance (all stable
    ones for the first 4 events, for the event displays), plus the outcome of every check."""
    import pickle
    def slim(rec):
        r = {k: v for k, v in rec.items() if k not in ("fulltruth", "gencands")}
        if "pf" in r: r["pf"] = {k: [x for x in v if abs(x[1]) < 2.5] for k, v in r["pf"].items()}
        return r
    H = {ev: [p for p in hep.get(ev, []) if p[1] == 1 and (i < 4 or (math.hypot(p[2], p[3]) > 0.5 and abs(kin(*p[2:])[1]) < 2.5))] for i, ev in enumerate(common)}
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(dict(gen=gen, results=RESULTS, common=common, atlas={ev: slim(atlas[ev]) for ev in common}, cms={ev: slim(cms[ev]) for ev in common}, hep=H), f, protocol=4)
    print(f"summary: {path}")

def merge(files, outdir):
    """Second step at scale: tally the checks of the per-seed passes and make the plots from their summaries."""
    import pickle
    nev, seeds, common, atlas, cms, hep, res = 0, [], [], {}, {}, {}, {}
    for f in files:
        with open(f, "rb") as fh: s = pickle.load(fh)
        nev += s["gen"]["nevents"]; seeds += s["gen"]["seeds"]; common += s["common"]; atlas.update(s["atlas"]); cms.update(s["cms"]); hep.update(s["hep"])
        for name, passed, detail in s["results"]: res.setdefault(name, []).append((s["gen"]["seeds"], passed, detail))
    print(f"merged {len(files)} summaries: seeds {min(seeds)}..{max(seeds)}, {nev} events generated, {len(common)} joined\n")
    ok = len(common) == nev
    for name, r in res.items():
        bad = [x for x in r if not x[1]]; ok &= not bad
        print(f"[{'PASS' if not bad else 'FAIL'}] {name}: {len(r) - len(bad)}/{len(r)} samples" + (f"; failed e.g. {[(b[0], b[2]) for b in bad[:3]]}" if bad else ""))
    os.makedirs(outdir, exist_ok=True)
    make_plots(sorted(common), atlas, cms, hep, outdir)
    print(f"\nplots: {outdir}/")
    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED"); return 0 if ok else 2

def make_plots(common, atlas, cms, hep, outdir):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    # (a) truth identity: stable truth e/mu pT ATLAS vs CMS; hard-process HepMC stable-particle summary
    tA, tC = [], []
    for ev in common:
        for p, s, pt, eta, phi in atlas[ev]["truth"]:
            if s != 1 or abs(p) not in (11, 13): continue
            m = [pt2 for q, s2, pt2, eta2, phi2 in cms[ev]["truth"] if q == p and s2 == 1 and dr((pt, eta, phi), (pt2, eta2, phi2)) < 0.01]
            if m: tA.append(pt); tC.append(m[0])
    # (b) reco vs reco: ΔR-matched electrons and jets between the two experiments (same event => same coordinates)
    eA, eC, jA, jC, nJ, met = [], [], [], [], [], []
    for ev in common:
        for a, c in match_objects(atlas[ev]["ele"], cms[ev]["ele"], 0.1): eA.append(a[0]); eC.append(c[0])
        for a, c in match_objects([j for j in atlas[ev]["jets"] if j[0] > 20], [j for j in cms[ev]["jets"] if j[0] > 20], 0.3): jA.append(a[0]); jC.append(c[0])
        nJ.append((sum(j[0] > 25 for j in atlas[ev]["jets"]), sum(j[0] > 25 for j in cms[ev]["jets"])))
        if atlas[ev]["met"] is not None: met.append((atlas[ev]["met"], cms[ev]["met"]))
    # (c) reco vs truth: electron response per experiment against HepMC stable electrons
    rA, rC = [], []
    for ev in common:
        te = [kin(px, py, pz, e) for pdg, st, px, py, pz, e in hep.get(ev, []) if st == 1 and abs(pdg) == 11 and math.hypot(px, py) > 5]
        for a, t in match_objects(atlas[ev]["ele"], te, 0.1): rA.append(a[0] / t[0])
        for c, t in match_objects(cms[ev]["ele"], te, 0.1): rC.append(c[0] / t[0])

    fig, ax = plt.subplots(2, 3, figsize=(15, 9))
    def ident(axis, x, y, title, unit="GeV"):
        if len(x) > 20000:   # a scatter saturates: log-density instead
            from matplotlib.colors import LogNorm
            lim = [0, 1.1 * max(max(x), max(y))]; axis.hist2d(x, y, bins=120, range=[lim, lim], norm=LogNorm(), cmap="Blues")
        else: axis.scatter(x, y, s=18)
        if x: lim = [0, 1.1 * max(max(x), max(y))]; axis.plot(lim, lim, "k--", lw=0.8); axis.set_xlim(lim); axis.set_ylim(lim)
        axis.set_xlabel(f"ATLAS [{unit}]"); axis.set_ylabel(f"CMS [{unit}]"); axis.set_title(f"{title} (n={len(x)})")
    ident(ax[0, 0], tA, tC, "truth e/μ pT (PHYSLITE vs GenPart)")
    ident(ax[0, 1], eA, eC, "reco electron pT (ΔR<0.1 matched)")
    ident(ax[0, 2], jA, jC, "reco jet pT >20 GeV (ΔR<0.3 matched)")
    ident(ax[1, 0], [m[0] for m in met], [m[1] for m in met], "MET (ATLAS core-term sum vs PuppiMET)")
    ax[1, 1].hist([[n[0] for n in nJ], [n[1] for n in nJ]], bins=range(0, 13), label=["ATLAS", "CMS"], alpha=0.7)
    ax[1, 1].set_xlabel("N jets (pT>25 GeV)"); ax[1, 1].set_ylabel("events"); ax[1, 1].legend(); ax[1, 1].set_title("jet multiplicity")
    ax[1, 2].hist([rA, rC], bins=np.linspace(0.7, 1.3, 25), label=[f"ATLAS (n={len(rA)})", f"CMS (n={len(rC)})"], alpha=0.7)
    ax[1, 2].set_xlabel("reco e pT / truth e pT"); ax[1, 2].set_ylabel("electrons"); ax[1, 2].legend(); ax[1, 2].set_title("electron response vs HepMC truth")
    fig.suptitle(f"{TITLE or os.path.basename(os.path.dirname(outdir))}: same {len(common)} Pythia8 events through ATLAS and CMS full simulation")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "compare.png"), dpi=110)

    # particle-flow level (ATLAS FlowElements from pflow.root vs CMS packed PF candidates from NanoAOD), same fiducial cuts on both:
    # charged: pT > 0.5 GeV, |eta| < 2.5 (tracker); neutral: pT > 1 GeV, |eta| < 2.5 — so that the two algorithms' different
    # low-pT/forward thresholds do not dominate the comparison.
    pfev = [ev for ev in common if "pf" in atlas[ev] and "pf" in cms[ev]]
    if pfev:
        sel_ch = lambda l: [p for p, e in l if p > 0.5 and abs(e) < 2.5]; sel_ne = lambda l: [p for p, e in l if p > 1.0 and abs(e) < 2.5]
        nch = [(len(sel_ch(atlas[ev]["pf"]["ch"])), len(sel_ch(cms[ev]["pf"]["ch"]))) for ev in pfev]
        sch = [(sum(sel_ch(atlas[ev]["pf"]["ch"])), sum(sel_ch(cms[ev]["pf"]["ch"]))) for ev in pfev]
        nne = [(len(sel_ne(atlas[ev]["pf"]["ne"])), len(sel_ne(cms[ev]["pf"]["ne"]))) for ev in pfev]
        sne = [(sum(sel_ne(atlas[ev]["pf"]["ne"])), sum(sel_ne(cms[ev]["pf"]["ne"]))) for ev in pfev]
        # truth reference: stable charged particles in the same acceptance
        tch = {ev: [kin(px, py, pz, e)[0] for pdg, st, px, py, pz, e in hep.get(ev, []) if st == 1 and abs(pdg) in (211, 321, 2212, 11, 13)
                    and kin(px, py, pz, e)[0] > 0.5 and abs(kin(px, py, pz, e)[1]) < 2.5] for ev in pfev}
        fig, ax = plt.subplots(2, 3, figsize=(15, 9))
        ident(ax[0, 0], [n[0] for n in nch], [n[1] for n in nch], "N charged PF (pT>0.5, |η|<2.5)", unit="count")
        ident(ax[0, 1], [n[0] for n in sch], [n[1] for n in sch], "Σ pT charged PF (pT>0.5, |η|<2.5)")
        ident(ax[0, 2], [n[0] for n in nne], [n[1] for n in nne], "N neutral PF (pT>1, |η|<2.5)", unit="count")
        # the decisive panel: charged + neutral together = the visible energy each experiment reconstructed
        tot = [(sum(p for p, e in atlas[ev]["pf"]["ch"] if abs(e) < 2.5) + sum(p for p, e in atlas[ev]["pf"]["ne"] if abs(e) < 2.5),
                sum(p for p, e in cms[ev]["pf"]["ch"] if abs(e) < 2.5) + sum(p for p, e in cms[ev]["pf"]["ne"] if abs(e) < 2.5)) for ev in pfev]
        ident(ax[1, 0], [t[0] for t in tot], [t[1] for t in tot], "Σ pT charged+neutral PF (|η|<2.5): same energy, different split")
        ms = 18 if len(pfev) < 5000 else 2
        ax[1, 1].scatter([len(tch[ev]) for ev in pfev], [n[0] for n in nch], label="ATLAS charged FE", s=ms)
        ax[1, 1].scatter([len(tch[ev]) for ev in pfev], [n[1] for n in nch], label="CMS charged PF", s=ms, marker="x")
        m = max([len(tch[ev]) for ev in pfev] + [1]); ax[1, 1].plot([0, m], [0, m], "k--", lw=0.8)
        ax[1, 1].set_xlabel("HepMC stable charged (pT>0.5, |η|<2.5)"); ax[1, 1].set_ylabel("reco charged"); ax[1, 1].legend(); ax[1, 1].set_title("charged multiplicity vs truth")
        bins = np.logspace(-0.3, 2.5, 30)
        ax[1, 2].hist([[p for ev in pfev for p in sel_ch(atlas[ev]["pf"]["ch"])], [p for ev in pfev for p in sel_ch(cms[ev]["pf"]["ch"])], [p for ev in pfev for p in tch[ev]]],
                      bins=bins, histtype="step", label=["ATLAS charged FE", "CMS charged PF", "HepMC stable charged"])
        ax[1, 2].set_xscale("log"); ax[1, 2].set_xlabel("charged pT [GeV]"); ax[1, 2].set_ylabel("candidates"); ax[1, 2].legend(); ax[1, 2].set_title("charged PF pT spectrum")
        fig.suptitle(f"{TITLE or os.path.basename(os.path.dirname(outdir))}: particle-flow level, {len(pfev)} events")
        fig.tight_layout(); fig.savefig(os.path.join(outdir, "compare_pflow.png"), dpi=110); plt.close(fig)

    make_marginals(common, atlas, cms, hep, outdir, plt)

    # per-event event display: truth stable particles vs reco jets/electrons for both experiments (first 4 events)
    for ev in common[:4]:
        fig, axs = plt.subplots(1, 2, figsize=(13, 5), sharex=True, sharey=True)
        st = [(kin(px, py, pz, e), pdg) for pdg, s, px, py, pz, e in hep.get(ev, []) if s == 1 and math.hypot(px, py) > 1 and abs(pdg) not in (12, 14, 16)]
        for axis, name, rec in [(axs[0], "ATLAS", atlas[ev]), (axs[1], "CMS", cms[ev])]:
            axis.scatter([k[1] for k, p in st], [k[2] for k, p in st], s=[3 * k[0] for k, p in st], c="lightgray", label="HepMC stable (size ∝ pT)")
            for j in rec["jets"]:
                if j[0] > 20: axis.add_patch(plt.Circle((j[1], j[2]), 0.4, fill=False, color="C0")); axis.text(j[1], j[2] + 0.45, f"{j[0]:.0f}", color="C0", fontsize=8, ha="center")
            for e in rec["ele"]: axis.plot(e[1], e[2], "r*", ms=12); axis.text(e[1], e[2] - 0.35, f"e {e[0]:.0f}", color="r", fontsize=8, ha="center")
            axis.set_xlim(-5, 5); axis.set_ylim(-math.pi, math.pi); axis.set_xlabel("η"); axis.set_ylabel("φ")
            axis.set_title(f"{name} reco, {ev[0]} event {ev[1]}: jets pT>20 GeV (circles), electrons (stars)")
        axs[0].legend(loc="upper left", fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(outdir, f"event_{ev[0]}_{ev[1]}.png"), dpi=100); plt.close(fig)

def has_partner(A, B, maxdr):
    """One flag per object of A: does it have a partner in B within maxdr (greedy in pT, each B used once)?"""
    flags, used = [False] * len(A), set()
    for i in sorted(range(len(A)), key=lambda k: -A[k][0]):
        best = min(((dr(A[i], y), j) for j, y in enumerate(B) if j not in used), default=None)
        if best and best[0] < maxdr: flags[i] = True; used.add(best[1])
    return flags

def make_marginals(common, atlas, cms, hep, outdir, plt):
    """The scatter panels of compare.png only contain objects found by BOTH experiments and ΔR-matched, so they cannot show what
    one detector reconstructs and the other does not. Here: every object of each experiment (marginals), the reconstruction
    efficiency of each against the common HepMC truth, and how often an object of one experiment has a partner in the other."""
    CA, CC = "C0", "C1"   # ATLAS, CMS: same colours in every panel; truth is always the grey fill
    ETA_E, PT_PROBE, PT_PARTNER = 2.5, 30., 15.
    eA, eC, jA, jC, tE, metA, metC, genmet, nE, nch = [], [], [], [], [], [], [], [], [], []
    for ev in common:
        a, c = atlas[ev], cms[ev]
        te = [k for k in (kin(px, py, pz, e) for pdg, st, px, py, pz, e in hep.get(ev, []) if st == 1 and abs(pdg) == 11) if k[0] > 5 and abs(k[1]) < ETA_E]
        tE += [(t, fa, fc) for t, fa, fc in zip(te, has_partner(te, a["ele"], 0.1), has_partner(te, c["ele"], 0.1))]
        eA += list(zip(a["ele"], has_partner(a["ele"], c["ele"], 0.1))); eC += list(zip(c["ele"], has_partner(c["ele"], a["ele"], 0.1)))
        # jets: probe above PT_PROBE, partner above PT_PARTNER, so that the partner's own pT threshold does not fake an inefficiency
        pa, pc = [j for j in a["jets"] if j[0] > PT_PROBE], [j for j in c["jets"] if j[0] > PT_PROBE]
        jA += list(zip(pa, has_partner(pa, [j for j in c["jets"] if j[0] > PT_PARTNER], 0.3)))
        jC += list(zip(pc, has_partner(pc, [j for j in a["jets"] if j[0] > PT_PARTNER], 0.3)))
        if a["met"] is not None: metA.append(a["met"])
        metC.append(c["met"]); nE.append((len(a["ele"]), len(c["ele"])))
        if c.get("genmet") is not None: genmet.append(c["genmet"])
        if "pf" in a and "pf" in c:
            sel = lambda l: sum(1 for p, e in l if p > 0.5 and abs(e) < 2.5)
            nch.append((sel(a["pf"]["ch"]), sel(c["pf"]["ch"]), sum(1 for pdg, st, px, py, pz, e in hep.get(ev, []) if st == 1 and abs(pdg) in (211, 321, 2212, 11, 13)
                                                                     and kin(px, py, pz, e)[0] > 0.5 and abs(kin(px, py, pz, e)[1]) < 2.5)))

    def marg(axis, series, bins, xlabel, title, truth=None, logy=False):
        if truth: axis.hist(truth[0], bins=bins, color="0.85", label=f"{truth[1]} (n={len(truth[0])})")
        for d, l, c in series: axis.hist(d, bins=bins, histtype="step", lw=1.6, color=c, label=f"{l} (n={len(d)})")
        if logy: axis.set_yscale("log")
        axis.set_xlabel(xlabel); axis.set_ylabel("entries"); axis.set_title(title); axis.legend(fontsize=8)
    def frac(axis, series, bins, xlabel, title, ylabel):
        bins = np.asarray(bins, dtype=float); ctr = np.sqrt(bins[:-1] * bins[1:]) if bins[0] > 0 else 0.5 * (bins[:-1] + bins[1:])
        for k, (pts, l, c, mk) in enumerate(series):   # pts = [(x, passed)]
            n, _ = np.histogram([x for x, p in pts], bins); kk, _ = np.histogram([x for x, p in pts if p], bins); ok = n > 0
            eff = kk[ok] / n[ok]; err = np.sqrt(np.clip(eff * (1 - eff), 0, None) / n[ok])
            tot = sum(p for x, p in pts) / max(len(pts), 1)
            axis.errorbar(ctr[ok], eff, yerr=err, fmt=mk + "-", ms=5, lw=1.2, color=c, label=f"{l}: {100 * tot:.1f}% of {len(pts)}")
        axis.set_ylim(0, 1.08); axis.axhline(1, color="0.7", lw=0.8); axis.grid(axis="y", color="0.9")
        if bins[0] > 0: axis.set_xscale("log"); axis.xaxis.set_minor_formatter(plt.NullFormatter())
        axis.set_xlabel(xlabel); axis.set_ylabel(ylabel); axis.set_title(title); axis.legend(fontsize=8, loc="lower right")

    fig, ax = plt.subplots(3, 4, figsize=(21, 14))
    ptb = np.logspace(math.log10(5), math.log10(500), 25); effb = [5, 7, 10, 15, 20, 30, 40, 60, 100, 500]; etab = np.linspace(-2.5, 2.5, 26)
    # electrons
    marg(ax[0, 0], [([e[0] for e, m in eA], "ATLAS", CA), ([e[0] for e, m in eC], "CMS", CC)], ptb, "electron pT [GeV]", "all reco electrons: pT",
         truth=([t[0] for t, fa, fc in tE], "HepMC e, pT>5, |η|<2.5"), logy=True); ax[0, 0].set_xscale("log")
    marg(ax[0, 1], [([e[1] for e, m in eA], "ATLAS", CA), ([e[1] for e, m in eC], "CMS", CC)], etab, "electron η", "all reco electrons: η",
         truth=([t[1] for t, fa, fc in tE], "HepMC e, pT>5"))
    frac(ax[0, 2], [([(t[0], fa) for t, fa, fc in tE], "ATLAS", CA, "o"), ([(t[0], fc) for t, fa, fc in tE], "CMS", CC, "s")], effb,
         "truth electron pT [GeV]", "electron efficiency vs the common truth (ΔR<0.1)", "fraction of truth electrons reconstructed")
    cat = [sum(fa and fc for t, fa, fc in tE), sum(fa and not fc for t, fa, fc in tE), sum(fc and not fa for t, fa, fc in tE), sum(not fa and not fc for t, fa, fc in tE)]
    bars = ax[0, 3].bar(["both", "ATLAS only", "CMS only", "neither"], cat, color=["0.45", CA, CC, "0.85"], width=0.6)
    for b, v in zip(bars, cat): ax[0, 3].text(b.get_x() + b.get_width() / 2, v, f"{v}\n{100 * v / max(sum(cat), 1):.1f}%", ha="center", va="bottom", fontsize=9)
    ax[0, 3].set_ylim(0, 1.2 * max(cat + [1])); ax[0, 3].set_ylabel("truth electrons"); ax[0, 3].set_title(f"who reconstructs each truth electron (n={sum(cat)})")
    # jets
    jb = np.logspace(math.log10(PT_PROBE), math.log10(1000), 25); jeb = [30, 40, 50, 70, 100, 150, 250, 1000]; jetab = np.linspace(-5, 5, 26)
    marg(ax[1, 0], [([j[0] for j, m in jA], "ATLAS", CA), ([j[0] for j, m in jC], "CMS", CC)], jb, "jet pT [GeV]", f"all reco jets, pT>{PT_PROBE:.0f} GeV: pT", logy=True); ax[1, 0].set_xscale("log")
    marg(ax[1, 1], [([j[1] for j, m in jA], "ATLAS", CA), ([j[1] for j, m in jC], "CMS", CC)], jetab, "jet η", f"all reco jets, pT>{PT_PROBE:.0f} GeV: η")
    part = f"fraction with a partner (ΔR<0.3, pT>{PT_PARTNER:.0f}) in the other experiment"
    frac(ax[1, 2], [([(j[0], m) for j, m in jA], "ATLAS jets found by CMS", CA, "o"), ([(j[0], m) for j, m in jC], "CMS jets found by ATLAS", CC, "s")], jeb, "jet pT [GeV]", "jets: cross-experiment match vs pT", part)
    frac(ax[1, 3], [([(j[1], m) for j, m in jA], "ATLAS jets found by CMS", CA, "o"), ([(j[1], m) for j, m in jC], "CMS jets found by ATLAS", CC, "s")], np.linspace(-5, 5, 11), "jet η", "jets: cross-experiment match vs η", part)
    # event level
    marg(ax[2, 0], [(metA, "ATLAS core-term sum", CA), (metC, "CMS PuppiMET", CC)], np.linspace(0, 100, 26), "MET [GeV]", "MET", truth=(genmet, "CMS GenMET") if genmet else None, logy=True)
    marg(ax[2, 1], [([n[0] for n in nE], "ATLAS", CA), ([n[1] for n in nE], "CMS", CC)], np.arange(-0.5, 7.5, 1), "reco electrons per event", "electron multiplicity")
    frac(ax[2, 2], [([(e[0], m) for e, m in eA], "ATLAS e found by CMS", CA, "o"), ([(e[0], m) for e, m in eC], "CMS e found by ATLAS", CC, "s")], effb, "reco electron pT [GeV]",
         "electrons: cross-experiment match vs pT", "fraction with a partner (ΔR<0.1) in the other experiment")
    if nch: marg(ax[2, 3], [([n[0] for n in nch], "ATLAS charged FE", CA), ([n[1] for n in nch], "CMS charged PF", CC)], np.linspace(0, 160, 41), "charged candidates per event (pT>0.5, |η|<2.5)",
                 "charged particle-flow multiplicity", truth=([n[2] for n in nch], "HepMC stable charged"))
    else: ax[2, 3].axis("off")
    fig.suptitle(f"{TITLE or os.path.basename(os.path.dirname(outdir))}: marginals and efficiencies, all objects (matched or not), {len(common)} events")
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(os.path.join(outdir, "compare_marginals.png"), dpi=100); plt.close(fig)

    pct = lambda pts: f"{100 * sum(m for o, m in pts) / max(len(pts), 1):.1f}% of {len(pts)}"
    print("\nunmatched objects (not visible in the scatter plots):")
    print(f"  truth electrons (pT>5, |η|<{ETA_E}): {sum(cat)}  -> both {cat[0]}, ATLAS only {cat[1]}, CMS only {cat[2]}, neither {cat[3]}")
    print(f"  reco electrons with a partner in the other experiment: ATLAS {pct(eA)}, CMS {pct(eC)}")
    print(f"  reco jets pT>{PT_PROBE:.0f} with a partner pT>{PT_PARTNER:.0f} in the other experiment: ATLAS {pct(jA)}, CMS {pct(jC)}")

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", help="one or more sample directories (each with gen/, atlas/, cms/); globs allowed")
    ap.add_argument("--out", help="plot directory (default: <dir>/compare for one dir, <parent>/compare_all for several)")
    ap.add_argument("--summary", help="also write a slim pickle of what the plots need (first step at scale)")
    ap.add_argument("--no-plots", action="store_true", help="checks (and --summary) only")
    ap.add_argument("--merge", action="store_true", help="the arguments are --summary files: tally their checks and make the plots (needs --out)")
    ap.add_argument("--title", help="plot title prefix (default: name of the directory above the plot directory)")
    args = ap.parse_args(); TITLE = args.title
    dirs = sorted(x for pat in args.dirs for x in (glob.glob(pat) or [pat]))
    if args.merge and not args.out: ap.error("--merge needs --out")
    sys.exit(merge(dirs, args.out) if args.merge else main(dirs, args.out, args.summary, not args.no_plots))
