#!/bin/bash
# Generate the same sample twice and assert the HepMC files are byte-identical.
set -eo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T="${TMPDIR:-/tmp}/gen_repro_$$"; mkdir -p "$T"
"${REPO}/gen/run_gen.sh" zee_13p6TeV 4242 50 "$T/a" >/dev/null
"${REPO}/gen/run_gen.sh" zee_13p6TeV 4242 50 "$T/b" >/dev/null
cmp "$T/a/events.hepmc3" "$T/b/events.hepmc3" && cmp "$T/a/events.hepmc2" "$T/b/events.hepmc2" \
  && echo "OK: generator is reproducible (identical HepMC3 and HepMC2 output for same seed)"
rm -rf "$T"
