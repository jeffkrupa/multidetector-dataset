#!/usr/bin/env python
"""Dump ATLAS particle-flow objects (and the full truth record) from an AOD into a flat ROOT TTree.

    python atlas/dump_pflow.py AOD.pool.root pflow.root        # run inside the Athena environment (asetup)

Why: PHYSLITE/PHYS do not keep FlowElements, and the AOD's static aux store is not uproot-readable.
Output tree "pflow", one entry per event, std::vector branches (GeV):
  event                                   event number (same as the HepMC event number)
  fe_ch_pt/eta/phi/e/charge               JetETMissChargedParticleFlowObjects (charged FlowElements: tracks matched to calo)
  fe_ne_pt/eta/phi/e                      JetETMissNeutralParticleFlowObjects (neutral FlowElements: calo clusters after subtraction)
  tp_pdg/status/pt/eta/phi/e              TruthParticles (full HepMC record as stored by ATLAS)
Reading the xAOD is done in C++ through ROOT's interpreter (xAOD::TEvent) because the PyROOT template call is fragile.
"""
import sys, ROOT
ROOT.xAOD.Init().ignore()
ROOT.gInterpreter.Declare(r'''
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODPFlow/FlowElementContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
int dump_pflow(const char* in, const char* out) {
  auto f = TFile::Open(in); if (!f) return 1;
  xAOD::TEvent ev(xAOD::TEvent::kClassAccess); if (ev.readFrom(f).isFailure()) return 2;
  TFile o(out, "RECREATE"); TTree t("pflow", "ATLAS particle-flow objects and truth (GeV)");
  ULong64_t event = 0; t.Branch("event", &event);
  std::vector<float> chpt, cheta, chphi, che, chq, nept, neeta, nephi, nee, tppt, tpeta, tpphi, tpe; std::vector<int> tppdg, tpst;
  t.Branch("fe_ch_pt",&chpt); t.Branch("fe_ch_eta",&cheta); t.Branch("fe_ch_phi",&chphi); t.Branch("fe_ch_e",&che); t.Branch("fe_ch_charge",&chq);
  t.Branch("fe_ne_pt",&nept); t.Branch("fe_ne_eta",&neeta); t.Branch("fe_ne_phi",&nephi); t.Branch("fe_ne_e",&nee);
  t.Branch("tp_pdg",&tppdg); t.Branch("tp_status",&tpst); t.Branch("tp_pt",&tppt); t.Branch("tp_eta",&tpeta); t.Branch("tp_phi",&tpphi); t.Branch("tp_e",&tpe);
  const Long64_t n = ev.getEntries();
  for (Long64_t i = 0; i < n; ++i) {
    ev.getEntry(i);
    const xAOD::EventInfo* ei = nullptr; ev.retrieve(ei, "EventInfo").ignore(); event = ei->eventNumber();
    const xAOD::FlowElementContainer *ch = nullptr, *ne = nullptr; const xAOD::TruthParticleContainer* tp = nullptr;
    ev.retrieve(ch, "JetETMissChargedParticleFlowObjects").ignore();
    ev.retrieve(ne, "JetETMissNeutralParticleFlowObjects").ignore();
    ev.retrieve(tp, "TruthParticles").ignore();
    chpt.clear(); cheta.clear(); chphi.clear(); che.clear(); chq.clear(); nept.clear(); neeta.clear(); nephi.clear(); nee.clear();
    tppdg.clear(); tpst.clear(); tppt.clear(); tpeta.clear(); tpphi.clear(); tpe.clear();
    for (auto p : *ch) { chpt.push_back(p->pt()/1e3); cheta.push_back(p->eta()); chphi.push_back(p->phi()); che.push_back(p->e()/1e3); chq.push_back(p->charge()); }
    for (auto p : *ne) { nept.push_back(p->pt()/1e3); neeta.push_back(p->eta()); nephi.push_back(p->phi()); nee.push_back(p->e()/1e3); }
    for (auto p : *tp) { tppdg.push_back(p->pdgId()); tpst.push_back(p->status()); tppt.push_back(p->pt()/1e3); tpeta.push_back(p->eta()); tpphi.push_back(p->phi()); tpe.push_back(p->e()/1e3); }
    t.Fill();
  }
  o.Write(); o.Close(); printf("dump_pflow: %lld events -> %s\n", n, out); return 0;
}''')
sys.exit(ROOT.dump_pflow(sys.argv[1], sys.argv[2]))
