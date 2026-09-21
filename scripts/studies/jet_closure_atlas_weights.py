import sys, glob, numpy as np, uproot, awkward as ak
sys.path.insert(0, sys.argv[3]); from recluster import antikt
d, N = sys.argv[1], int(sys.argv[2]); sa = lambda c, v: f"{c}Aux./{c}Aux.{v}"; ch, ne, trk = "JetETMissChargedParticleFlowObjects", "JetETMissNeutralParticleFlowObjects", "InDetTrackParticles"
t = uproot.open(glob.glob(d + "/atlas/DAOD_PHYSLITE.*.root")[0])["CollectionTree"]; KEYS = set(t.keys())
_sa = sa; sa = lambda c, v: _sa(c, v) if _sa(c, v) in KEYS else f"{c}AuxDyn.{v}"   # static or dynamic aux variable
A = t.arrays([sa(ch, v) for v in ("pt", "eta", "phi", "m", "IsInDenseEnvironment", "TracksExpectedEnergyDeposit", "chargedObjectLinks")] + [sa(ne, v) for v in ("pt", "eta", "phi", "m", "CENTER_MAG")] + ["PrimaryVerticesAuxDyn.x", "PrimaryVerticesAuxDyn.y"] + [sa(trk, v) for v in ("z0", "vz", "theta")] +
             ["PrimaryVerticesAuxDyn.z", "PrimaryVerticesAuxDyn.vertexType"] + [f"AnalysisJetsAuxDyn.{v}" for v in ("pt", "eta", "phi", "JetConstitScaleMomentum_pt")], entry_stop=N)
res = {k: [] for k in ("raw list", "weights + vertex", "weights + vertex + neutral origin corr.")}; stat = dict(n=0, dense=0, w0=0, failv=0, ptw=0., pt=0.)
for i in range(N):
    g = lambda c, v: np.asarray(A[sa(c, v)][i], float); pt, eta, phi, m = (g(ch, v) for v in ("pt", "eta", "phi", "m")); dense = g(ch, "IsInDenseEnvironment") > 0; eexp = g(ch, "TracksExpectedEnergyDeposit")
    e = np.sqrt((pt * np.cosh(eta))**2 + m**2); eop = eexp / e; w = np.where(pt < 30e3, 1., np.where(pt < 60e3, eop + (1 - (pt - 30e3) / 30e3) * (1 - eop), eop)); w = np.where(dense, w - eop, w); w = np.where(pt > 100e3, 0., w); w = np.clip(w, 0, None)
    idx = np.asarray(ak.fill_none(ak.firsts(A[sa(ch, "chargedObjectLinks")][i]["m_persIndex"]), -1)); pvt = np.asarray(A["PrimaryVerticesAuxDyn.vertexType"][i]); pvz = np.asarray(A["PrimaryVerticesAuxDyn.z"][i])[pvt == 1]
    if len(pvz) and np.all(idx >= 0):
        z0, vz, th = (np.asarray(A[sa(trk, v)][i], float)[idx] for v in ("z0", "vz", "theta")); okv = np.abs((z0 + vz - pvz[0]) * np.sin(th)) < 2.0
    else: okv = np.ones(len(pt), bool)
    stat["n"] += len(pt); stat["dense"] += dense.sum(); stat["w0"] += (w < 1e-6).sum(); stat["failv"] += (~okv).sum(); stat["pt"] += pt.sum(); stat["ptw"] += (pt * w * okv).sum()
    npt, neta, nphi, nm, mag = (g(ne, v) for v in ("pt", "eta", "phi", "m", "CENTER_MAG")); k0 = npt > 0
    # origin correction: keep the cluster energy, point it from the primary vertex instead of the detector centre
    pv = np.array([np.asarray(A["PrimaryVerticesAuxDyn." + c][i])[pvt == 1][0] for c in "xyz"]) if len(pvz) else np.zeros(3)
    pos = mag[:, None] * np.stack([np.cos(nphi) / np.cosh(neta), np.sin(nphi) / np.cosh(neta), np.tanh(neta)], axis=1) - pv[None, :]; rr = np.linalg.norm(pos, axis=1)
    ceta = np.where(mag > 0, np.arctanh(np.clip(pos[:, 2] / np.maximum(rr, 1e-9), -0.999999, 0.999999)), neta); cphi = np.where(mag > 0, np.arctan2(pos[:, 1], pos[:, 0]), nphi)
    en = np.sqrt((npt * np.cosh(neta))**2 + nm**2); cpt = np.sqrt(np.maximum(en**2 - nm**2, 0)) / np.cosh(ceta)
    for name, cw, (upt, ueta, uphi) in (("raw list", np.ones(len(pt)), (npt, neta, nphi)), ("weights + vertex", w * okv, (npt, neta, nphi)), ("weights + vertex + neutral origin corr.", w * okv, (cpt, ceta, cphi))):
        k = cw > 1e-6; J = antikt(np.concatenate([pt[k] * cw[k], upt[k0]]) / 1e3, np.concatenate([eta[k], ueta[k0]]), np.concatenate([phi[k], uphi[k0]]), np.concatenate([m[k] * cw[k], np.maximum(nm[k0], 0)]) / 1e3)
        for jpt, jraw, jeta, jphi in zip(A["AnalysisJetsAuxDyn.pt"][i] / 1e3, A["AnalysisJetsAuxDyn.JetConstitScaleMomentum_pt"][i] / 1e3, A["AnalysisJetsAuxDyn.eta"][i], A["AnalysisJetsAuxDyn.phi"][i]):
            if jpt < 30 or abs(jeta) > 2 or len(J[0]) == 0: continue
            dr = np.hypot(J[1] - jeta, (J[2] - jphi + np.pi) % (2 * np.pi) - np.pi); q = int(np.argmin(dr)); res[name].append(J[0][q] / jraw if dr[q] < 0.1 else np.nan)
print(f"charged flow objects: {stat['n']}; in dense environment {100*stat['dense']/stat['n']:.1f}%; weight 0 {100*stat['w0']/stat['n']:.1f}%; fail |z0 sin(theta)|<2mm {100*stat['failv']/stat['n']:.1f}%; charged pT kept after weights+vertex {100*stat['ptw']/stat['pt']:.1f}%")
print(f"{'recipe applied to the stored list':42s} {'jets':>5s} {'found':>7s} {'median':>8s} {'16%':>7s} {'84%':>7s} {'within 1%':>10s} {'within 5%':>10s}")
for name, r in res.items():
    r = np.array(r); ok = ~np.isnan(r); q = np.percentile(r[ok], [16, 50, 84]); print(f"{name:42s} {len(r):5d} {100*ok.mean():6.1f}% {q[1]:8.3f} {q[0]:7.3f} {q[2]:7.3f} {100*np.mean(np.abs(r[ok]-1)<0.01):9.1f}% {100*np.mean(np.abs(r[ok]-1)<0.05):9.1f}%")
