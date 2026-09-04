"""Neutrino Discovery Platform (NDP) — theory -> generator -> surrogate detector -> data.

The package is organised along the pipeline:

    theory/      a theorist's model spec, realised as a truth-level event sample
                 (GENIE run, reweighted reference MC, external sample) or directly as
                 a cross-section vector (a shipped generator curve)
    adapters/    readers that turn experiment / generator files into the common
                 `events.TruthTable` (MINERvA MasterAnaDev AnaTuples, GENIE gst)
    channels/    channel manifests: signal definition, phase space, observables,
                 binning, selection, data references, normalisation constants
    surrogate/   detector surrogates p(reco | truth): a binned response
                 (efficiency x migration) and a parametric smearing model
    compare/     the two comparison modes — unfolded space (published d2sigma +
                 covariance) and folded space (surrogate-smeared prediction vs
                 reconstructed data counts)
    pipeline.py  orchestration: one call produces a manifest-backed run directory

Units inside the platform are GeV (energies, momenta), GeV^2 (Q^2), mm (vertices),
cm^2 (cross sections). Adapters convert on the way in.
"""
__version__ = "0.1.0"
