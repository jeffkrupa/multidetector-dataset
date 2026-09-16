import FWCore.ParameterSet.Config as cms
def customise(process):
    """GEN,SIM from an external HepMC2 ASCII file (cmsDriver was given a dummy fragment).
    Replaces the source and turns 'generator' into a relabeler that presents the
    standard generator products (+ empty GenRunInfoProduct / GenLumiInfoHeader)."""
    process.source = cms.Source("MCFileSource",
        fileNames = cms.untracked.vstring("file:events.hepmc"),
        firstLuminosityBlockForEachRun = cms.untracked.VLuminosityBlockID())
    # cmsDriver already inserted process.generator at the front of every path;
    # reassigning the label swaps the dummy Pythia8 EDFilter for our producer.
    process.generator = cms.EDProducer("HepMCSourceRelabeler",
        src = cms.InputTag("source", "generator"))
    return process
