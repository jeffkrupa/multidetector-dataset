# Athena Gen_tf.py job configuration: read externally generated HepMC (2 or 3) ASCII events into EVNT.
# Directory name (999002) is the local DSID; file name must be mc.<Generator>_<...>.py (Generator short name "HepMC").
# The HepMC file is picked up by Gen_tf.py --inputGeneratorFile=<name>.events from the run directory and
# symlinked to events.hepmc, which TruthIO/HepMCReadFromFile reads (HepMC2/3 auto-detected via HepMC3::deduce_reader).
evgenConfig.description = "Externally generated Pythia8 HepMC events (shared ATLAS/CMS sample)"
evgenConfig.keywords = ["SM"]
evgenConfig.contact = ["jeffkrupa@gmail.com"]
evgenConfig.nEventsPerJob = 100000     # upper bound; actual count is --maxEvents
evgenConfig.inputFilesPerJob = 1
evgenConfig.tune = "Monash 2013 (external Pythia8)"
include("TruthIO/HepMCReadFromFile_Common.py")

# ATLAS must take every event CMS takes. TestHepMC drops ~1 event in 10^4 (a decay vertex > 1 m from the origin, or
# > 100 mm transverse; Pythia energy-momentum imbalance of tens of GeV); Gen_tf then reads past the end of the HepMC
# file to make up the count and aborts (~13% of 1000-event jobs). Switch off only those rejections, keep the rest
# (NaN, tachyons, unknown PDG IDs, unstable particles without a decay vertex, ...).
if hasattr(testSeq, "TestHepMC"):
    testSeq.TestHepMC.VtxDisplacedTest = False
    testSeq.TestHepMC.EnergyImbalanceTest = False
    testSeq.TestHepMC.MomImbalanceTest = False
