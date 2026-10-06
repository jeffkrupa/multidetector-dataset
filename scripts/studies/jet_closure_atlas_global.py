import sys, glob, numpy as np, uproot
sys.path.insert(0, sys.argv[2]); from recluster import antikt
d = sys.argv[1]; N = 100
tf = uproot.open(d + "/atlas/DAOD_FTAG1.out.pool.root")["CollectionTree"]; K = set(tf.keys())
def br(c, v):   # FTAG1 writes the Global containers as AuxDyn branches
    for k in (f"{c}AuxDyn.{v}", f"{c}Aux./{c}Aux.{v}"):
        if k in K: return k
    raise KeyError(c + "." + v)
GC, GN = "GlobalChargedParticleFlowObjects", "GlobalNeutralParticleFlowObjects"
jetvars = {v: f"AntiKt4EMPFlowJetsAuxDyn.{v}" for v in ("pt", "eta", "phi", "JetConstitScaleMomentum_pt")}
have_jets = all(j in K for j in jetvars.values())
print("Global containers present:", all(br(c, "pt") for c in (GC, GN)), "| FTAG1 jet collections:", sorted({k.split("Aux")[0] for k in K if "Jets" in k and "Aux" in k})[:6])
F = tf.arrays([br(c, v) for c in (GC, GN) for v in ("pt", "eta", "phi", "m")] + [br(GC, "IsInDenseEnvironment"), br(GC, "TracksExpectedEnergyDeposit"), "EventInfoAuxDyn.eventNumber"] + (list(jetvars.values()) if have_jets else []), entry_stop=N, library="np")
tp = uproot.open(d + "/atlas/DAOD_PHYSLITE.out.pool.root")["CollectionTree"]
P = tp.arrays([f"AnalysisJetsAuxDyn.{v}" for v in ("pt", "eta", "phi", "JetConstitScaleMomentum_pt")] + [f"AnalysisElectronsAuxDyn.{v}" for v in ("pt", "eta", "phi")] + [f"AnalysisMuonsAuxDyn.{v}" for v in ("pt", "eta", "phi")] + [f"AnalysisPhotonsAuxDyn.{v}" for v in ("pt", "eta", "phi")] + ["EventInfoAuxDyn.eventNumber"], entry_stop=N, library="np")
ip = {int(e): i for i, e in enumerate(P["EventInfoAuxDyn.eventNumber"])}
res = {k: [] for k in ("Global, as stored", "Global + weights", "Global + weights + e/mu/gamma added back")}; ng = nj = 0
for i in range(N):
    j = ip[int(F["EventInfoAuxDyn.eventNumber"][i])]
    g = lambda c, v: np.asarray(F[br(c, v)][i], float); pt, eta, phi, m = (g(GC, v) for v in ("pt", "eta", "phi", "m")); npt, neta, nphi, nm = (g(GN, v) for v in ("pt", "eta", "phi", "m"))
    dense = g(GC, "IsInDenseEnvironment") > 0; eop = g(GC, "TracksExpectedEnergyDeposit") / np.sqrt((pt * np.cosh(eta))**2 + m**2)
    w = np.where(pt < 30e3, 1., np.where(pt < 60e3, eop + (1 - (pt - 30e3) / 30e3) * (1 - eop), eop)); w = np.clip(np.where(dense, w - eop, w), 0, None); w = np.where(pt > 100e3, 0, w)
    lep = [np.concatenate([np.asarray(P[f"Analysis{c}AuxDyn.{v}"][j], float) for c in ("Electrons", "Muons", "Photons")]) for v in ("pt", "eta", "phi")]
    k0 = npt > 0; ng += len(pt) + k0.sum(); nj += 0
    for name, cw, extra in (("Global, as stored", np.ones(len(pt)), None), ("Global + weights", w, None), ("Global + weights + e/mu/gamma added back", w, lep)):
        k = cw > 1e-6; parts = [np.concatenate([pt[k] * cw[k], npt[k0]]), np.concatenate([eta[k], neta[k0]]), np.concatenate([phi[k], nphi[k0]]), np.concatenate([m[k] * cw[k], np.maximum(nm[k0], 0)])]
        if extra is not None: parts = [np.concatenate([parts[0], extra[0]]), np.concatenate([parts[1], extra[1]]), np.concatenate([parts[2], extra[2]]), np.concatenate([parts[3], np.zeros(len(extra[0]))])]
        J = antikt(parts[0] / 1e3, parts[1], parts[2], parts[3] / 1e3)
        for jpt, jraw, jeta, jphi in zip(P["AnalysisJetsAuxDyn.pt"][j] / 1e3, P["AnalysisJetsAuxDyn.JetConstitScaleMomentum_pt"][j] / 1e3, P["AnalysisJetsAuxDyn.eta"][j], P["AnalysisJetsAuxDyn.phi"][j]):
            if jpt < 30 or abs(jeta) > 2 or len(J[0]) == 0: continue
            dr = np.hypot(J[1] - jeta, (J[2] - jphi + np.pi) % (2 * np.pi) - np.pi); q = int(np.argmin(dr)); res[name].append(J[0][q] / jraw if dr[q] < 0.1 else np.nan)
print(f"{'Global particles -> anti-kT R=0.4, vs PHYSLITE AnalysisJets (pT>30, |eta|<2)':62s} {'jets':>5s} {'found':>7s} {'median':>8s} {'16%':>7s} {'84%':>7s} {'within 1%':>10s} {'within 5%':>10s}")
for name, r in res.items():
    r = np.array(r); ok = ~np.isnan(r); q = np.percentile(r[ok], [16, 50, 84]); print(f"{name:62s} {len(r):5d} {100*ok.mean():6.1f}% {q[1]:8.3f} {q[0]:7.3f} {q[2]:7.3f} {100*np.mean(np.abs(r[ok]-1)<0.01):9.1f}% {100*np.mean(np.abs(r[ok]-1)<0.05):9.1f}%")
# how much is removed by overlap removal: compare object counts with the JetETMiss containers of the PHYSLITE file
Q = tp.arrays(["JetETMissChargedParticleFlowObjectsAux./JetETMissChargedParticleFlowObjectsAux.pt", "JetETMissNeutralParticleFlowObjectsAux./JetETMissNeutralParticleFlowObjectsAux.pt"], entry_stop=N, library="np")
nje = sum(len(x) for x in Q["JetETMissChargedParticleFlowObjectsAux./JetETMissChargedParticleFlowObjectsAux.pt"]) + sum((np.asarray(x) > 0).sum() for x in Q["JetETMissNeutralParticleFlowObjectsAux./JetETMissNeutralParticleFlowObjectsAux.pt"])
print(f"objects per event: JetETMiss {nje/N:.0f}, Global {ng/N:.0f}  ({100*(1-ng/nje):.1f}% removed by overlap removal)")
