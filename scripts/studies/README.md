One-off studies behind docs/RECO_DIFFERENCES.md, kept as run (LCG python; on a non-EL9 host run inside the EL9 image, see
scripts/compare.sh). Not part of the workflow.

  jet_closure_cms_atlas_raw.py <seed dir> <N>        anti-kT on the stored particles vs the stored jets, both experiments
  jet_closure_atlas_weights.py <seed dir> <N> <dir of jet_closure_cms_atlas_raw.py renamed recluster.py>
                                                     same with the ATLAS WeightPFOTool weights and vertex association
  electron_id_scan.py <prod dir> <process> ...       electron ID working points of both experiments vs prompt / non-prompt truth
  electron_eff_vs_pt_isolation.py <prod dir> <process> ...   prompt-electron efficiency at fixed pT and truth isolation
