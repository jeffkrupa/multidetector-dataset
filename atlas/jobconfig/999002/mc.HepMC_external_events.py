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
