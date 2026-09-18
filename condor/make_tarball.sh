#!/bin/bash
# Pack the repository (sources only: no output/, no build areas, no .git) into condor/repo.tar.gz for shipping with each
# job, and record provenance (commit, dirty files) inside it.
set -eo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "${REPO}"
{ echo "commit: $(git rev-parse HEAD 2>/dev/null || echo none)"; echo "date: $(date -u +%FT%TZ)"; echo "uncommitted:"; git status --short 2>/dev/null; } > PROVENANCE
tar czf condor/repo.tar.gz --exclude=output --exclude=cms/build --exclude=gen/build --exclude=.git --exclude=condor/repo.tar.gz --exclude='*.log' .
echo "condor/repo.tar.gz: $(du -h condor/repo.tar.gz | cut -f1), $(head -1 PROVENANCE)"
