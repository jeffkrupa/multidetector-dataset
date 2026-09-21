import sys, glob, numpy as np, uproot
def antikt(pt, eta, phi, m, R=0.4, ptmin=10.):
    """Plain N^3 anti-kT with E-scheme recombination (fine for ~200 particles)."""
    px, py = pt * np.cos(phi), pt * np.sin(phi); pz = pt * np.sinh(eta); E = np.sqrt(px**2 + py**2 + pz**2 + m**2)
    P = [np.array(v) for v in zip(px, py, pz, E)]; jets = []
    while P:
        A = np.array(P); ptj = np.hypot(A[:, 0], A[:, 1]); y = 0.5 * np.log((A[:, 3] + A[:, 2]) / np.maximum(A[:, 3] - A[:, 2], 1e-12)); ph = np.arctan2(A[:, 1], A[:, 0])
        diB = 1. / ptj**2; n = len(P)
        if n > 1:
            dphi = (ph[:, None] - ph[None, :] + np.pi) % (2 * np.pi) - np.pi; d = np.minimum(diB[:, None], diB[None, :]) * ((y[:, None] - y[None, :])**2 + dphi**2) / R**2
            d[np.arange(n), np.arange(n)] = np.inf; i, j = np.unravel_index(np.argmin(d), d.shape); dmin = d[i, j]
        else: dmin = np.inf
        k = int(np.argmin(diB))
        if diB[k] <= dmin: jets.append(P.pop(k))
        else:
            new = P[i] + P[j]
            for q in sorted((i, j), reverse=True): P.pop(q)
            P.append(new)
    J = np.array([j for j in jets if np.hypot(j[0], j[1]) > ptmin]).reshape(-1, 4)
    p = np.sqrt(J[:, 0]**2 + J[:, 1]**2 + J[:, 2]**2); return np.hypot(J[:, 0], J[:, 1]), np.arctanh(J[:, 2] / np.maximum(p, 1e-9)), np.arctan2(J[:, 1], J[:, 0])
def compare(tag, stored, reco, out):
    for spt, sraw, seta, sphi in stored:
        if spt < 30 or abs(seta) > 2.0 or len(reco[0]) == 0: continue
        dr = np.hypot(reco[1] - seta, (reco[2] - sphi + np.pi) % (2 * np.pi) - np.pi); k = int(np.argmin(dr))
        out.append((dr[k] < 0.1, reco[0][k] / sraw if dr[k] < 0.1 else np.nan, reco[0][k] / spt if dr[k] < 0.1 else np.nan))
d = sys.argv[1]; N = int(sys.argv[2]); sa = lambda c, v: f"{c}Aux./{c}Aux.{v}"; ch, ne = "JetETMissChargedParticleFlowObjects", "JetETMissNeutralParticleFlowObjects"
A = uproot.open(glob.glob(d + "/atlas/DAOD_PHYSLITE.*.root")[0])["CollectionTree"].arrays([sa(c, v) for c in (ch, ne) for v in ("pt", "eta", "phi", "m")] +
      [f"AnalysisJetsAuxDyn.{v}" for v in ("pt", "eta", "phi", "JetConstitScaleMomentum_pt")], entry_stop=N, library="np")
C = uproot.open(d + "/cms/step4_nano.root")["Events"].arrays([f"PFCands_{v}" for v in ("pt", "eta", "phi", "mass", "puppiWeight")] + [f"Jet_{v}" for v in ("pt", "eta", "phi", "rawFactor")], entry_stop=N, library="np")
ra, rc, rc0 = [], [], []
for i in range(N):
    f = lambda c, v: np.asarray(A[sa(c, v)][i], float); cat = lambda v: np.concatenate([f(ch, v), f(ne, v)])
    compare("A", zip(A["AnalysisJetsAuxDyn.pt"][i] / 1e3, A["AnalysisJetsAuxDyn.JetConstitScaleMomentum_pt"][i] / 1e3, A["AnalysisJetsAuxDyn.eta"][i], A["AnalysisJetsAuxDyn.phi"][i]),
            antikt(cat("pt") / 1e3, cat("eta"), cat("phi"), np.maximum(cat("m"), 0) / 1e3), ra)
    w = np.asarray(C["PFCands_puppiWeight"][i], float); g = lambda v: np.asarray(C[f"PFCands_{v}"][i], float); st = list(zip(C["Jet_pt"][i], C["Jet_pt"][i] * (1 - C["Jet_rawFactor"][i]), C["Jet_eta"][i], C["Jet_phi"][i]))
    k = w > 0; compare("C", st, antikt(g("pt")[k] * w[k], g("eta")[k], g("phi")[k], g("mass")[k] * w[k]), rc); compare("C0", st, antikt(g("pt"), g("eta"), g("phi"), g("mass")), rc0)
for lab, r in (("ATLAS: stored particles -> anti-kT R=0.4, vs AnalysisJets", ra), ("CMS: PFCands x puppiWeight -> anti-kT R=0.4, vs Jet", rc), ("CMS: PFCands without puppi weights, vs Jet", rc0)):
    r = np.array(r, float); m = r[:, 0] > 0; q = lambda x: np.nanpercentile(x, [16, 50, 84])
    print(f"{lab}\n   stored jets pT>30, |eta|<2: {len(r)}; reclustered jet within dR<0.1: {100*m.mean():.1f}%")
    print("   reclustered pT / stored constituent-scale (raw) pT: median %.3f  [16%%, 84%%] = [%.3f, %.3f];  within 1%%: %.1f%%" % (q(r[m, 1])[1], q(r[m, 1])[0], q(r[m, 1])[2], 100 * np.mean(np.abs(r[m, 1] - 1) < 0.01)))
    print("   reclustered pT / stored calibrated pT:              median %.3f" % q(r[m, 2])[1])
