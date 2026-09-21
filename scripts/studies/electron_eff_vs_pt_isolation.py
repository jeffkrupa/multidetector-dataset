import sys, glob, math, uproot, numpy as np
def dr(e1, p1, e2, p2): return math.hypot(e1 - e2, (p1 - p2 + math.pi) % (2 * math.pi) - math.pi)
PTB = [10, 20, 30, 45, 70, 1e9]; SELS = ["ATLAS all", "ATLAS LHMedium", "CMS all", "CMS mvaNoIso_WP90", "CMS cutBased>=2"]
res = {}
for proc in sys.argv[2:]:
    R = res[proc] = {}; comp = {"prompt": 0, "nonprompt": 0}
    for d in sorted(glob.glob(f"{sys.argv[1]}/{proc}/seed*"))[:8]:
        p = "AnalysisElectronsAuxDyn."
        A = uproot.open(glob.glob(d + "/atlas/DAOD_PHYSLITE.*.root")[0])["CollectionTree"].arrays([p + v for v in ["pt", "eta", "phi", "DFCommonElectronsLHMedium"]] + ["EventInfoAuxDyn.eventNumber"], library="np")
        C = uproot.open(d + "/cms/step4_nano.root")["Events"].arrays(["event", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_cutBased", "Electron_mvaNoIso_WP90",
              "GenPart_pdgId", "GenPart_status", "GenPart_statusFlags", "GenPart_eta", "GenPart_phi", "GenCands_pdgId", "GenCands_pt", "GenCands_eta", "GenCands_phi"], library="np")
        ia = {int(e): i for i, e in enumerate(A["EventInfoAuxDyn.eventNumber"])}
        for j, ev in enumerate(C["event"]):
            i = ia[int(ev)]
            pr = [(e, f) for g, s, fl, e, f in zip(C["GenPart_pdgId"][j], C["GenPart_status"][j], C["GenPart_statusFlags"][j], C["GenPart_eta"][j], C["GenPart_phi"][j]) if abs(g) == 11 and s == 1 and (fl & 1 or fl & 4)]
            gc = list(zip(C["GenCands_pdgId"][j], C["GenCands_pt"][j], C["GenCands_eta"][j], C["GenCands_phi"][j]))
            sets = {"ATLAS all": [(e, f) for pt, e, f in zip(A[p + "pt"][i], A[p + "eta"][i], A[p + "phi"][i])],
                    "ATLAS LHMedium": [(e, f) for e, f, m in zip(A[p + "eta"][i], A[p + "phi"][i], A[p + "DFCommonElectronsLHMedium"][i]) if m],
                    "CMS all": list(zip(C["Electron_eta"][j], C["Electron_phi"][j])),
                    "CMS mvaNoIso_WP90": [(e, f) for e, f, w in zip(C["Electron_eta"][j], C["Electron_phi"][j], C["Electron_mvaNoIso_WP90"][j]) if w],
                    "CMS cutBased>=2": [(e, f) for e, f, w in zip(C["Electron_eta"][j], C["Electron_phi"][j], C["Electron_cutBased"][j]) if w >= 2]}
            for g, pt, e, f in gc:
                if abs(g) != 11 or pt < 10 or abs(e) > 2.5: continue
                isp = any(dr(e, f, e2, f2) < 0.01 for e2, f2 in pr); comp["prompt" if isp else "nonprompt"] += 1
                if not isp: continue
                iso = sum(q for g2, q, e2, f2 in gc if abs(g2) not in (12, 14, 16) and 1e-4 < dr(e, f, e2, f2) < 0.3) / pt   # truth relative isolation, same for both detectors
                b = next(k for k in range(5) if pt < PTB[k + 1]); env = "isolated (truth relIso<0.1)" if iso < 0.1 else "busy (truth relIso>0.1)"
                for key in (("pt", b), ("env", env), ("envpt", env, b >= 2)):
                    for s in SELS:
                        c = R.setdefault((key, s), [0, 0]); c[1] += 1; c[0] += any(dr(e, f, e2, f2) < 0.1 for e2, f2 in sets[s])
    print(f"\n== {proc}: truth electrons pT>10: {comp['prompt']} prompt, {comp['nonprompt']} non-prompt ({100*comp['nonprompt']/sum(comp.values()):.0f}% non-prompt)")
hdr = f"{'':42s}" + "".join(f"{s:>20s}" for s in SELS)
for title, keys in [("PROMPT-electron efficiency vs truth pT", [("pt", b) for b in range(5)]), ("... vs truth isolation (all pT>10)", [("env", e) for e in ("isolated (truth relIso<0.1)", "busy (truth relIso>0.1)")]),
                    ("... vs truth isolation, pT>30 only", [("envpt", e, True) for e in ("isolated (truth relIso<0.1)", "busy (truth relIso>0.1)")])]:
    print("\n" + title); print(hdr)
    for k in keys:
        for proc in sys.argv[2:]:
            R = res[proc]; n = R.get((k, SELS[0]), [0, 0])[1]
            lab = (f"pT {PTB[k[1]]}-{PTB[k[1]+1] if PTB[k[1]+1] < 1e8 else 'inf'}" if k[0] == "pt" else k[1]) + f"  {proc.split('_')[0]} (n={n})"
            print(f"  {lab:40s}" + "".join(f"{(100*R[(k,s)][0]/R[(k,s)][1] if (k,s) in R and R[(k,s)][1] else float('nan')):19.1f}%" for s in SELS))
