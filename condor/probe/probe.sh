#!/bin/bash
# Worker-node probe: what does an HTCondor node at CERN actually see?
echo "host=$(hostname) date=$(date) user=$(id -un) cpus=${_CONDOR_JOB_CPUS:-?} scratch=${_CONDOR_SCRATCH_DIR:-?}"
echo "OS: $(grep PRETTY /etc/os-release)"
echo "== kerberos"; klist 2>&1 | head -3
echo "== /eos mounted?"; mount | grep -c " /eos"; ls /eos 2>&1 | head -3; ls /eos/home-j/jeffkrup/agentic/multidetector-dataset 2>&1 | head -5
echo "== cvmfs"; ls /cvmfs/cms.cern.ch/common/cmssw-el9 /cvmfs/atlas.cern.ch/repo/containers/fs/singularity/x86_64-almalinux9 2>&1 | head -3
echo "== apptainer"; which apptainer; apptainer --version 2>&1
echo "== xrootd to EOS with kerberos"; xrdfs eosuser.cern.ch ls /eos/home-j/jeffkrup/agentic/multidetector-dataset 2>&1 | head -5
echo "== xrdcp test"; echo probe > probe_$(hostname).txt && xrdcp -f probe_$(hostname).txt root://eosuser.cern.ch//eos/home-j/jeffkrup/agentic/multidetector-dataset/output/condor/probe_$(hostname).txt 2>&1 | tail -1 && echo xrdcp_ok
echo "== apptainer exec test"; apptainer exec -B /cvmfs /cvmfs/atlas.cern.ch/repo/containers/fs/singularity/x86_64-almalinux9 cat /etc/os-release 2>&1 | grep PRETTY
echo "== disk"; df -h . | tail -1
echo "probe done"
