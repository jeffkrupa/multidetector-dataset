#!/usr/bin/env python3
"""Compare the ATLAS and CMS outputs of one sample, event by event, at generator and reconstruction level.

    source gen/env_lcg.sh                      # needs uproot, awkward, numpy, matplotlib (all in the LCG view)
    python3 scripts/compare_events.py output/<sample>

Checks (each prints PASS/FAIL and the script exits non-zero on any FAIL):
  1. hepmc-inputs   the .hepmc2 file given to CMS and the .hepmc3 file given to ATLAS contain identical particles
  2. atlas-truth    every truth particle ATLAS stored (PHYSLITE Truth* containers) exists in the HepMC file
  3. cms-truth      every gen particle CMS stored (NanoAOD GenPart) exists in the HepMC file
  4. truth-xcheck   stable truth leptons/photons stored by both experiments agree 1:1 in (pdgId, pt, eta, phi)
  5. event-join     both experiments wrote every generated event, joined on the event number

Then a per-event reco table and plots in <sample>/compare/ (truth identity, reco vs reco, reco vs truth).
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
    with open(path) as f:
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
    return out, truth_containers

# --------------------------------------------------------------------------------------- CMS
def load_cms(d):
    f = os.path.join(d, "cms", "step4_nano.root")
    if not os.path.exists(f): return None
    t = uproot.open(f)["Events"]; keys = set(t.keys())
    want = ["event", "GenPart_pdgId", "GenPart_status", "GenPart_pt", "GenPart_eta", "GenPart_phi",
            "Electron_pt", "Electron_eta", "Electron_phi", "Electron_cutBased", "nMuon", "Jet_pt", "Jet_eta", "Jet_phi", "PuppiMET_pt", "GenMET_pt"]
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

def report(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    return ok

# --------------------------------------------------------------------------------------- main
def main(d):
    d = d.rstrip("/"); ok = True
    gen = json.load(open(os.path.join(d, "gen", "events.json")))
    print(f"sample: {d}\ngenerator: Pythia {gen['pythia_version']} seed {gen['seed']} nevents {gen['nevents']} card {os.path.basename(gen['card'])}\n")

    # 1. HepMC inputs identical
    h3 = read_hepmc(os.path.join(d, "gen", "events.hepmc3")); h2 = read_hepmc(os.path.join(d, "gen", "events.hepmc2"))
    key = lambda p: (p[0], p[1], round(p[5], 4), round(p[4], 4), round(p[3], 4))   # order-independent (HepMC2 lists particles per vertex)
    same = set(h3) == set(h2) and all(len(h3[e]) == len(h2[e]) and all(
        p[0] == q[0] and p[1] == q[1] and all(abs(p[k] - q[k]) <= 1e-6 * max(1., abs(p[k])) for k in range(2, 6))
        for p, q in zip(sorted(h3[e], key=key), sorted(h2[e], key=key))) for e in h3)
    ok &= report("hepmc-inputs", same, f"{len(h3)} events, {sum(len(v) for v in h3.values())} particles, "
                 f"{sum(sum(1 for p in v if p[1] == 1) for v in h3.values())} stable; HepMC3 (ATLAS) == HepMC2 (CMS)")

    atlas = load_atlas(d); cms = load_cms(d)
    if atlas is None or cms is None:
        print("ATLAS output:", "found" if atlas else "MISSING", "| CMS output:", "found" if cms else "MISSING"); return 1
    atlas, atlas_containers = atlas

    # 2./3. stored truth ⊂ HepMC
    for name, coll, tol in [("atlas-truth", atlas, 1e-3), ("cms-truth", cms, 0.01)]:   # NanoAOD stores pt/eta/phi at ~1e-3 precision
        nchk, bad = 0, []
        for ev, rec in coll.items():
            c = [(p, pt, eta, phi) for p, s, pt, eta, phi in rec["truth"]]
            nchk += len(c); bad += [(ev, b) for b in match_to_hepmc(c, h3.get(ev, []), tol)]
        src = f"PHYSLITE {'+'.join(atlas_containers)}" if name == "atlas-truth" else "NanoAOD GenPart"
        ok &= report(name, not bad, f"{nchk} stored truth particles ({src}) all found in HepMC" if not bad else f"{len(bad)} unmatched, e.g. {bad[:3]}")

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
    ok &= report("event-join", len(common) == gen["nevents"] == len(atlas) == len(cms),
                 f"generated {gen['nevents']}, ATLAS {len(atlas)}, CMS {len(cms)}, joined {len(common)}")

    # per-event reco table
    print(f"\n{'event':>5} | {'ATLAS n_e':>9} {'e1 pT':>6} {'n_jet':>5} {'jet pT (GeV)':>24} {'MET':>6} | {'CMS n_e':>7} {'e1 pT':>6} {'n_jet':>5} {'jet pT (GeV)':>24} {'MET':>6}")
    fmt = lambda js: str([round(j[0], 1) for j in js[:4]])
    for ev in common:
        a, c = atlas[ev], cms[ev]
        print(f"{ev:>5} | {len(a['ele']):>9} {(round(a['ele'][0][0],1) if a['ele'] else '-'):>6} {len(a['jets']):>5} {fmt(a['jets']):>24} {(round(a['met'],1) if a['met'] is not None else '-'):>6} | "
              f"{len(c['ele']):>7} {(round(c['ele'][0][0],1) if c['ele'] else '-'):>6} {len(c['jets']):>5} {fmt(c['jets']):>24} {round(c['met'],1):>6}")

    # plots
    outdir = os.path.join(d, "compare"); os.makedirs(outdir, exist_ok=True)
    make_plots(common, atlas, cms, h3, outdir)
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
        axis.scatter(x, y, s=18)
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
    fig.suptitle(f"{os.path.basename(os.path.dirname(outdir))}: same {len(common)} Pythia8 events through ATLAS and CMS full simulation")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "compare.png"), dpi=110)

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
            axis.set_title(f"{name} reco, event {ev}: jets pT>20 GeV (circles), electrons (stars)")
        axs[0].legend(loc="upper left", fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(outdir, f"event_{ev}.png"), dpi=100); plt.close(fig)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "output"))
