"""Derivation_tf.py postInclude: add constituent-level content to the DAOD_PHYSLITE output stream, so one ATLAS file
carries analysis objects (PHYSLITE) plus particle-flow objects before overlap removal, the full truth record and,
optionally, tracks and topoclusters. Everything is written the derivation way (dynamic aux branches), so uproot reads it.

    Derivation_tf.py --CA --formats PHYSLITE --postInclude "default:MultiDetPHYSLITE.AddMultiDetContent" ...

Branch naming in the output: the added containers keep their native aux store, split into per-variable branches
    <Container>Aux./<Container>Aux.<variable>      e.g. TruthParticlesAux./TruthParticlesAux.pdgId
(uproot reads these directly), while PHYSLITE's own slimmed containers use <Container>AuxDyn.<variable>.
When "tracks" is requested, PHYSLITE's InDetTrackParticles thinning tools are removed from the kernel so that every
reconstructed track is kept (PHYSLITE alone keeps only tracks attached to analysis objects).

Content groups (env ATLAS_EXTRA_CONTENT, default "pflow truth tracks clusters"):
  pflow    JetETMissCharged/NeutralParticleFlowObjects  (all FlowElements BEFORE e/gamma/mu/tau overlap removal; CMS-PF analogue)
  truth    TruthParticles, TruthVertices                (full generator record as stored by ATLAS)
  tracks   InDetTrackParticles                          (perigee parameters, fit quality, hit counts, truth link)
  clusters CaloCalTopoClusters                          (calibrated and raw 4-vectors, moments)
"""
import os
from AthenaCommon.Logging import logging
log = logging.getLogger("MultiDetPHYSLITE")

CONTENT = {
    "pflow": [
        "xAOD::FlowElementContainer#JetETMissChargedParticleFlowObjects",
        "xAOD::FlowElementAuxContainer#JetETMissChargedParticleFlowObjectsAux.",
        "xAOD::FlowElementContainer#JetETMissNeutralParticleFlowObjects",
        "xAOD::FlowElementAuxContainer#JetETMissNeutralParticleFlowObjectsAux.",
    ],
    "truth": [
        "xAOD::TruthParticleContainer#TruthParticles",
        "xAOD::TruthParticleAuxContainer#TruthParticlesAux.",
        "xAOD::TruthVertexContainer#TruthVertices",
        "xAOD::TruthVertexAuxContainer#TruthVerticesAux.",
    ],
    "tracks": [
        "xAOD::TrackParticleContainer#InDetTrackParticles",
        "xAOD::TrackParticleAuxContainer#InDetTrackParticlesAux.d0.z0.phi.theta.qOverP.vz.chiSquared.numberDoF"
        ".numberOfInnermostPixelLayerHits.numberOfPixelHits.numberOfSCTHits.numberOfTRTHits.numberOfPixelHoles.numberOfSCTHoles"
        ".truthParticleLink.truthMatchProbability.truthType.truthOrigin",
    ],
    "clusters": [
        "xAOD::CaloClusterContainer#CaloCalTopoClusters",
        "xAOD::CaloClusterAuxContainer#CaloCalTopoClustersAux.calE.calEta.calPhi.calM.rawE.rawEta.rawPhi.rawM.time.clusterSize"
        ".CENTER_LAMBDA.SECOND_R.SECOND_LAMBDA.EM_PROBABILITY.ENG_FRAC_EM.ISOLATION",
    ],
}

def AddMultiDetContent(flags, cfg):
    groups = os.environ.get("ATLAS_EXTRA_CONTENT", "pflow truth tracks clusters").split()
    items = [i for g in groups for i in CONTENT[g]]
    stream = None
    for alg in cfg.getEventAlgos():
        if alg.getType() == "AthenaOutputStream" and "PHYSLITE" in alg.getName():
            stream = alg; break
    if stream is None:
        names = [f"{a.getName()}({a.getType()})" for a in cfg.getEventAlgos() if a.getType() == "AthenaOutputStream"]
        raise RuntimeError(f"MultiDetPHYSLITE: PHYSLITE output stream not found; output streams present: {names}")
    new = [i for i in items if i not in stream.ItemList]
    stream.ItemList += new
    log.info("MultiDetPHYSLITE: added %d items (%s) to %s", len(new), " ".join(groups), stream.getName())
    if "tracks" in groups:
        # PHYSLITE thins InDetTrackParticles with several tools (tracks of electrons, muons, taus, di-taus, generic
        # object-associated tracks); a track survives if ANY tool keeps it, so all of them must go to keep every track.
        # egamma calo-cluster thinning (egammaClusters) is unrelated and stays.
        TRACK_THINNERS = ("TrackParticleThinning", "TPThinning", "TauJetThinning", "TauJets_MuonRMThinning", "DiTau", "MuonTP")
        removed = []
        for alg in cfg.getEventAlgos():
            if alg.getName() != "PHYSLITEKernel": continue
            keep = []
            for tool in alg.ThinningTools:
                if any(k in tool.getName() for k in TRACK_THINNERS): removed.append(tool.getName())
                else: keep.append(tool)
            alg.ThinningTools = keep
            log.info("MultiDetPHYSLITE: thinning tools kept: %s", [t.getName() for t in keep])
        log.info("MultiDetPHYSLITE: removed track thinning tools so all InDetTrackParticles are kept: %s", removed)
