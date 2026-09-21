#!/usr/bin/env python3
"""Particle-level marginals: every particle-flow candidate of each experiment against the common truth, matched or not.

    scripts/compare.sh --particles --hist h/seed1000.npz <dir>/seed1000          # one pass per seed: histograms only
    scripts/compare.sh --particles --merge --out <plotdir> --title <name> h/*.npz  # add them up, compare_particles.png

ATLAS: JetETMiss{Charged,Neutral}ParticleFlowObjects of the DAOD_PHYSLITE (all flow elements, before overlap removal).
CMS:   PFCands of the NanoAOD (all packed PF candidates; pdgId 211/11/13 charged, 22 photon, 130 neutral hadron, 1/2 HF).
Truth: the stable generator particles, read from the CMS GenCands (compare_events.py verifies per seed that they are the
       HepMC stable particles with |eta|<5.9 and that ATLAS holds the same record). Neutrinos are excluded.
       With ParticleDecays:tau0Max = 10 mm, K0S, Lambda and the charged hyperons are stable here and decay in Geant4.
Matching is only used for the efficiency panels: a truth particle is "found" if a candidate of the right kind lies within
dR < 0.03 (charged) or 0.1 (photons); a charged candidate is "unmatched" if no charged truth particle (pT > 0.3) does.
"""
import argparse, glob, math, os, sys
import numpy as np
import uproot

CHARGED = (211, 321, 2212, 11, 13, 3112, 3222, 3312, 3334); NEUTRINOS = (12, 14, 16)
PT = np.logspace(math.log10(0.5), math.log10(200), 41); ETA = np.linspace(-5, 5, 51); ETAC = np.linspace(-2.5, 2.5, 26)
EPT = np.array([0.5, 0.7, 1, 1.5, 2, 3, 5, 10, 20, 50, 200.]); GPT = np.array([1, 1.5, 2, 3, 5, 10, 20, 50, 200.]); NCH = np.arange(-0.5, 200.5, 4); NNE = np.arange(-0.5, 80.5, 2)
BINS = dict(ch_pt=PT, ne_pt=PT, ch_eta=ETAC, ne_eta=ETA, flow_ch=ETA, flow_ne=ETA, n_ch=NCH, n_ne=NNE, eff_pt=EPT, eff_eta=ETAC, fake_pt=EPT, geff_pt=GPT)

def near(e1, p1, e2, p2, R):
    """For each (e1, p1): is there any (e2, p2) within dR < R?"""
    if len(e1) == 0 or len(e2) == 0: return np.zeros(len(e1), bool)
    dphi = (p1[:, None] - p2[None, :] + np.pi) % (2 * np.pi) - np.pi
    return (np.hypot(e1[:, None] - e2[None, :], dphi) < R).any(axis=1)

def one_seed(d):
    H = {}
    def add(name, x, w=None):
        kind = name.split(":")[1]; h = np.histogram(x, BINS[kind.replace("_num", "").replace("_den", "")], weights=w)[0]
        H[name] = H.get(name, 0) + h
    sa = lambda c, v: f"{c}Aux./{c}Aux.{v}"; ch, ne = "JetETMissChargedParticleFlowObjects", "JetETMissNeutralParticleFlowObjects"
    A = uproot.open(glob.glob(os.path.join(d, "atlas", "DAOD_PHYSLITE.*.pool.root"))[0])["CollectionTree"].arrays(
        [sa(c, v) for c in (ch, ne) for v in ("pt", "eta", "phi")] + ["EventInfoAuxDyn.eventNumber"], library="np")
    C = uproot.open(os.path.join(d, "cms", "step4_nano.root"))["Events"].arrays(
        ["event"] + [f"PFCands_{v}" for v in ("pt", "eta", "phi", "pdgId")] + [f"GenCands_{v}" for v in ("pt", "eta", "phi", "pdgId")], library="np")
    ia = {int(e): i for i, e in enumerate(A["EventInfoAuxDyn.eventNumber"])}; nev = 0
    for j, ev in enumerate(C["event"]):
        i = ia.get(int(ev))
        if i is None: continue
        nev += 1
        g, gpt, geta, gphi = (np.asarray(C[f"GenCands_{v}"][j]) for v in ("pdgId", "pt", "eta", "phi")); ag = np.abs(g)
        tch = np.isin(ag, CHARGED); tne = ~tch & ~np.isin(ag, NEUTRINOS); tga = ag == 22
        pid = np.abs(np.asarray(C["PFCands_pdgId"][j])); cch = np.isin(pid, (211, 11, 13))
        cand = {"ATLAS": {k: tuple(np.asarray(A[sa(c, v)][i], dtype=float) / (1000. if v == "pt" else 1.) for v in ("pt", "eta", "phi")) for k, c in (("ch", ch), ("ne", ne))},
                "CMS": {k: tuple(np.asarray(C[f"PFCands_{v}"][j], dtype=float)[m] for v in ("pt", "eta", "phi")) for k, m in (("ch", cch), ("ne", ~cch))},
                "truth": {k: (gpt[m], geta[m], gphi[m]) for k, m in (("ch", tch), ("ne", tne))}}
        for k, m in (("gamma", pid == 22), ("h0", pid == 130), ("hf", np.isin(pid, (1, 2)))):   # CMS neutral classes, for the eta panel
            add("CMS:ne_eta:" + k, np.asarray(C["PFCands_eta"][j])[m & (np.asarray(C["PFCands_pt"][j]) > 1)])
        for det, c in cand.items():
            (cpt, ceta, cphi), (npt, neta, nphi) = c["ch"], c["ne"]; cin = (np.abs(ceta) < 2.5) & (cpt > 0.5); nin = np.abs(neta) < 2.5
            add(det + ":ch_pt", cpt[cin]); add(det + ":ch_eta", ceta[cin]); add(det + ":ne_pt", npt[nin]); add(det + ":ne_eta", neta[npt > 1])
            add(det + ":flow_ch", ceta, cpt); add(det + ":flow_ne", neta, npt); add(det + ":n_ch", [cin.sum()]); add(det + ":n_ne", [(nin & (npt > 1)).sum()])
        (tpt, teta, tphi) = cand["truth"]["ch"]; tin = (np.abs(teta) < 2.5) & (tpt > 0.5); tloose = tpt > 0.3
        gin = tga & (np.abs(geta) < 2.5) & (gpt > 1)
        for det in ("ATLAS", "CMS"):
            (cpt, ceta, cphi), (npt, neta, nphi) = cand[det]["ch"], cand[det]["ne"]; cin = (np.abs(ceta) < 2.5) & (cpt > 0.5)
            found = near(teta[tin], tphi[tin], ceta, cphi, 0.03)
            add(det + ":eff_pt_den", tpt[tin]); add(det + ":eff_pt_num", tpt[tin][found]); hi = tpt[tin] > 1
            add(det + ":eff_eta_den", teta[tin][hi]); add(det + ":eff_eta_num", teta[tin][hi & found])
            real = near(ceta[cin], cphi[cin], teta[tloose], tphi[tloose], 0.03)
            add(det + ":fake_pt_den", cpt[cin]); add(det + ":fake_pt_num", cpt[cin][~real])
            gf = near(geta[gin], gphi[gin], neta, nphi, 0.1)
            add(det + ":geff_pt_den", gpt[gin]); add(det + ":geff_pt_num", gpt[gin][gf])
    H["nev"] = np.array([nev]); return H

def plot(H, outdir, title):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    nev = float(H["nev"][0]); COL = {"ATLAS": "C0", "CMS": "C1"}
    def marg(ax, kind, xlabel, ttl, tlabel, logx=False, logy=False, per_event=False, density_eta=False, extra=()):
        b = BINS[kind]; w = np.diff(b) if density_eta else 1.; s = 1. / nev if per_event else 1.
        ax.stairs(H["truth:" + kind] * s / w, b, fill=True, color="0.85", label=tlabel)
        for det in ("ATLAS", "CMS"): ax.stairs(H[det + ":" + kind] * s / w, b, color=COL[det], lw=1.6, label=det)
        for key, lab, ls in extra: ax.stairs(H[key] * s / w, b, color=COL["CMS"], lw=1, ls=ls, label=lab)
        if logx: ax.set_xscale("log")
        if logy: ax.set_yscale("log")
        ax.set_xlabel(xlabel); ax.set_title(ttl, fontsize=10); ax.legend(fontsize=8)
    def ratio(ax, kind, xlabel, ttl, ylabel, logx=True, ylim=(0, 1.08)):
        b = BINS[kind]; ctr = np.sqrt(b[:-1] * b[1:]) if logx else 0.5 * (b[:-1] + b[1:])
        for det, mk in (("ATLAS", "o"), ("CMS", "s")):
            n, k = H[f"{det}:{kind}_den"].astype(float), H[f"{det}:{kind}_num"].astype(float); ok = n > 0; e = k[ok] / n[ok]
            ax.errorbar(ctr[ok], e, yerr=np.sqrt(np.clip(e * (1 - e), 0, None) / n[ok]), fmt=mk + "-", ms=4, lw=1.2, color=COL[det], label=f"{det}: {100 * k.sum() / max(n.sum(), 1):.1f}% of {int(n.sum())}")
        if logx: ax.set_xscale("log"); ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.set_ylim(*ylim); ax.axhline(1, color="0.7", lw=0.8); ax.grid(axis="y", color="0.9"); ax.set_xlabel(xlabel); ax.set_ylabel(ylabel, fontsize=9); ax.set_title(ttl, fontsize=10); ax.legend(fontsize=8)
    fig, ax = plt.subplots(3, 4, figsize=(21, 14))
    marg(ax[0, 0], "ch_pt", "pT [GeV]", "charged candidates: pT  (|η|<2.5)", "truth stable charged", logx=True, logy=True); ax[0, 0].set_ylabel("candidates")
    marg(ax[0, 1], "ch_eta", "η", "charged candidates: η  (pT>0.5 GeV)", "truth stable charged"); ax[0, 1].set_ylabel("candidates")
    ratio(ax[0, 2], "eff_pt", "truth pT [GeV]", "charged: efficiency vs truth pT  (|η|<2.5, ΔR<0.03)", "fraction of truth charged particles found")
    ratio(ax[0, 3], "eff_eta", "truth η", "charged: efficiency vs truth η  (pT>1 GeV)", "fraction of truth charged particles found", logx=False)
    marg(ax[1, 0], "ne_pt", "pT [GeV]", "neutral candidates: pT  (|η|<2.5)", "truth stable neutral (γ, n, K0L, K0S, Λ)", logx=True, logy=True); ax[1, 0].set_ylabel("candidates")
    marg(ax[1, 1], "ne_eta", "η", "neutral candidates: η  (pT>1 GeV)", "truth stable neutral",
         extra=(("CMS:ne_eta:gamma", "CMS photons", "--"), ("CMS:ne_eta:h0", "CMS neutral hadrons", ":"), ("CMS:ne_eta:hf", "CMS HF", "-."))); ax[1, 1].set_ylabel("candidates")
    ratio(ax[1, 2], "geff_pt", "truth photon pT [GeV]", "photons: a neutral candidate within ΔR<0.1  (|η|<2.5)", "fraction of truth photons")
    ratio(ax[1, 3], "fake_pt", "candidate pT [GeV]", "charged candidates with no truth charged particle within ΔR<0.03", "fraction of candidates", ylim=(0, 0.5))
    marg(ax[2, 0], "flow_ch", "η", "charged energy flow: dΣpT/dη per event", "truth", per_event=True, density_eta=True); ax[2, 0].set_ylabel("GeV per unit η")
    marg(ax[2, 1], "flow_ne", "η", "neutral energy flow: dΣpT/dη per event", "truth", per_event=True, density_eta=True); ax[2, 1].set_ylabel("GeV per unit η")
    marg(ax[2, 2], "n_ch", "charged candidates per event (pT>0.5, |η|<2.5)", "charged multiplicity", "truth"); ax[2, 2].set_ylabel("events")
    marg(ax[2, 3], "n_ne", "neutral candidates per event (pT>1, |η|<2.5)", "neutral multiplicity", "truth"); ax[2, 3].set_ylabel("events")
    fig.suptitle(f"{title}: particle-level marginals, every particle-flow candidate (matched or not), {int(nev)} events")
    fig.tight_layout(rect=[0, 0, 1, 0.97]); os.makedirs(outdir, exist_ok=True); fig.savefig(os.path.join(outdir, "compare_particles.png"), dpi=100)
    print(f"{title}: {int(nev)} events -> {outdir}/compare_particles.png")
    for kind, lab in (("eff_pt", "truth charged (pT>0.5, |η|<2.5) found"), ("geff_pt", "truth photons (pT>1, |η|<2.5) with a neutral candidate"), ("fake_pt", "charged candidates without a truth particle")):
        print(f"  {lab}: " + ", ".join(f"{det} {100 * H[f'{det}:{kind}_num'].sum() / max(H[f'{det}:{kind}_den'].sum(), 1):.1f}%" for det in ("ATLAS", "CMS")))
    for kind, lab in (("flow_ch", "charged ΣpT per event, all η"), ("flow_ne", "neutral ΣpT per event, all η")):
        print(f"  {lab}: " + ", ".join(f"{det} {H[f'{det}:{kind}'].sum() / nev:.1f} GeV" for det in ("truth", "ATLAS", "CMS")))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="sample directories (seed*/), or .npz files with --merge")
    ap.add_argument("--hist", help="write the histograms of the given sample directories to this .npz")
    ap.add_argument("--merge", action="store_true", help="the arguments are --hist files: add them up and plot (needs --out)")
    ap.add_argument("--out", help="plot directory"); ap.add_argument("--title", default="")
    a = ap.parse_args(); paths = sorted(x for p in a.paths for x in (glob.glob(p) or [p])); H = {}
    for p in paths:
        h = dict(np.load(p)) if a.merge else one_seed(p.rstrip("/"))
        for k, v in h.items(): H[k] = H.get(k, 0) + v
    if a.hist: os.makedirs(os.path.dirname(os.path.abspath(a.hist)), exist_ok=True); np.savez_compressed(a.hist, **H); print(f"histograms: {a.hist} ({int(H['nev'][0])} events)")
    if a.out: plot(H, a.out, a.title or os.path.basename(os.path.dirname(os.path.abspath(paths[0]))))
