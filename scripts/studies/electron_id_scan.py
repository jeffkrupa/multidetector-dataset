import sys, glob, math, uproot, numpy as np
def dr(e1, p1, e2, p2): return math.hypot(e1 - e2, (p1 - p2 + math.pi) % (2 * math.pi) - math.pi)
ISO = "ptvarcone30_Nonprompt_All_MaxWeightTTVALooseCone_pt1000"
for proc in sys.argv[2:]:
    seeds = sorted(glob.glob(f"{sys.argv[1]}/{proc}/seed*"))[:4]
    selA = ["all AnalysisElectrons", "LHLoose", "LHMedium", "LHMedium + trkiso<0.15"]; selC = ["all NanoAOD Electrons", "cutBased>=1 (veto)", "cutBased>=2 (loose) [used now]", "mvaNoIso_WP90", "mvaNoIso_WP90 + relIso<0.15"]
    T = {"prompt": 0, "nonprompt": 0}; hit = {s: {"prompt": 0, "nonprompt": 0, "fake": 0, "n": 0} for s in selA + selC}; nev = 0
    for d in seeds:
        ta = uproot.open(glob.glob(d + "/atlas/DAOD_PHYSLITE.*.root")[0])["CollectionTree"]; p = "AnalysisElectronsAuxDyn."
        A = ta.arrays([p + v for v in ["pt", "eta", "phi", "DFCommonElectronsLHLoose", "DFCommonElectronsLHMedium", ISO]] + ["EventInfoAuxDyn.eventNumber"], library="np")
        C = uproot.open(d + "/cms/step4_nano.root")["Events"].arrays(["event", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_cutBased", "Electron_mvaNoIso_WP90", "Electron_pfRelIso03_all",
              "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags", "GenPart_pt", "GenPart_eta", "GenPart_phi", "GenCands_pdgId", "GenCands_pt", "GenCands_eta", "GenCands_phi"], library="np")
        ia = {int(e): i for i, e in enumerate(A["EventInfoAuxDyn.eventNumber"])}
        for j, ev in enumerate(C["event"]):
            i = ia.get(int(ev));
            if i is None: continue
            nev += 1
            pr = [(e, f) for g, s, fl, pt, e, f in zip(C["GenPart_pdgId"][j], C["GenPart_status"][j], C["GenPart_statusFlags"][j], C["GenPart_pt"][j], C["GenPart_eta"][j], C["GenPart_phi"][j]) if abs(g) == 11 and s == 1 and (fl & 1 or fl & 4)]
            truth = [(pt, e, f, "prompt" if any(dr(e, f, e2, f2) < 0.01 for e2, f2 in pr) else "nonprompt") for g, pt, e, f in zip(C["GenCands_pdgId"][j], C["GenCands_pt"][j], C["GenCands_eta"][j], C["GenCands_phi"][j]) if abs(g) == 11 and pt > 10 and abs(e) < 2.5]
            for t in truth: T[t[3]] += 1
            aa = [(pt / 1e3, e, f, bool(lo), bool(me), iso / pt) for pt, e, f, lo, me, iso in zip(*(A[p + v][i] for v in ["pt", "eta", "phi", "DFCommonElectronsLHLoose", "DFCommonElectronsLHMedium", ISO])) if pt > 1e4 and abs(e) < 2.5]
            cc = [(pt, e, f, cb, bool(w), iso) for pt, e, f, cb, w, iso in zip(*(C["Electron_" + v][j] for v in ["pt", "eta", "phi", "cutBased", "mvaNoIso_WP90", "pfRelIso03_all"])) if pt > 10 and abs(e) < 2.5]
            sets = {selA[0]: aa, selA[1]: [x for x in aa if x[3]], selA[2]: [x for x in aa if x[4]], selA[3]: [x for x in aa if x[4] and x[5] < 0.15],
                    selC[0]: cc, selC[1]: [x for x in cc if x[3] >= 1], selC[2]: [x for x in cc if x[3] >= 2], selC[3]: [x for x in cc if x[4]], selC[4]: [x for x in cc if x[4] and x[5] < 0.15]}
            for s, L in sets.items():
                hit[s]["n"] += len(L)
                for t in truth:
                    if any(dr(t[1], t[2], x[1], x[2]) < 0.1 for x in L): hit[s][t[3]] += 1
                hit[s]["fake"] += sum(1 for x in L if not any(dr(t[1], t[2], x[1], x[2]) < 0.1 for t in truth))
    print(f"\n== {proc}: {nev} events; truth electrons pT>10 |eta|<2.5: {T['prompt']} prompt (incl. from tau), {T['nonprompt']} non-prompt; reco electrons pT>10")
    print(f"   {'selection':34s} {'n reco':>7s} {'eff prompt':>11s} {'eff non-prompt':>15s} {'no truth e within 0.1':>22s}")
    for s in selA + selC:
        h = hit[s]; print(f"   {('ATLAS ' if s in selA else 'CMS   ') + s:34s} {h['n']:7d} {100*h['prompt']/max(T['prompt'],1):10.1f}% {100*h['nonprompt']/max(T['nonprompt'],1):14.1f}% {h['fake']:14d} ({1000*h['fake']/nev:.0f}/1k ev)")
