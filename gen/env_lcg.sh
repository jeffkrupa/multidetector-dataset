# Pinned LCG software stack for the generator step (Pythia8 + HepMC3 + gcc14, EL9).
# Source this; it is used both on a bare EL9 host with CVMFS and inside the EL9 apptainer image.
export LCG_RELEASE=LCG_110
export LCG_PLATFORM=x86_64-el9-gcc14-opt
export LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/${LCG_RELEASE}/${LCG_PLATFORM}
source "${LCG_VIEW}/setup.sh"
