import FWCore.ParameterSet.Config as cms

def allPF(process):
    """NanoAOD + ALL packed PF candidates (not only jet constituents), via the in-release PFNano (custom_btv_cff).
    Adds PFCands_* (pt, eta, phi, mass, pdgId, charge, puppiWeight, ...) for every packedPFCandidate, the
    jet<->candidate index tables (JetPFCands_*, FatJetPFCands_*) and the gen-level GenCands_* for all packedGenParticles."""
    from PhysicsTools.NanoAOD.custom_btv_cff import btvNano_switch, BTVCustomNanoAOD
    btvNano_switch.btvNano_addallPF_switch = True
    return BTVCustomNanoAOD(process)
