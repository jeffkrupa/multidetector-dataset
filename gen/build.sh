#!/bin/bash
# Build gen_hepmc against the pinned LCG view. Output: gen/build/gen_hepmc
# (paths are kept relative after cd: `mkdir -p` with absolute paths trips over EOS FUSE ACLs)
set -eo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source ./env_lcg.sh
mkdir -p build
PY8CFG="$(command -v pythia8-config)"
g++ -O2 -std=c++17 -o build/gen_hepmc gen_hepmc.cc \
    $(${PY8CFG} --cxxflags) $(${PY8CFG} --libs) \
    -I"${LCG_VIEW}/include" -L"${LCG_VIEW}/lib" -L"${LCG_VIEW}/lib64" -lHepMC3 \
    -Wl,-rpath,"${LCG_VIEW}/lib" -Wl,-rpath,"${LCG_VIEW}/lib64"
echo "built $(pwd)/build/gen_hepmc (Pythia8 $(${PY8CFG} --version 2>/dev/null || echo ?), ${LCG_RELEASE}/${LCG_PLATFORM})"
