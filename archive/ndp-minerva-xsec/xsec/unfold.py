"""D'Agostini iterative unfolding (RooUnfold kBayes equivalent).

Dimension-agnostic port of the frozen exploration-repo `dagostini_unfold`
(dsigma_dpt.py:334), operating on flat slot vectors with a count-conserving
migration oriented M[reco_slot, true_slot] (xsec.binning.migration_slots).

The response is the column-normalization of the migration,
R[reco, true] = M[reco, true] / Σ_reco M = P(reco | true); its true-column
sums are the efficiency numerator. Because the migration conserves counts
(under/overflow included), folding the true distribution through R reproduces
the reco distribution exactly, so MC self-closure is exact to machine
precision (see tests).

Error propagation uses the final-iteration Bayes matrix U applied to the
diagonal data covariance — the standard approximation (ignores the
inter-iteration feedback term), and is labelled as such by the caller.
"""
import numpy as np


def response_matrix(migration):
    """(response R[reco, true], colsum[true]) from migration M[reco, true].

    R columns are normalized to 1 where the true slot has events; empty true
    columns stay all-zero.
    """
    migration = np.asarray(migration, dtype=np.float64)
    colsum = migration.sum(axis=0)                       # eff numerator per true slot
    response = np.divide(migration, colsum,
                         out=np.zeros_like(migration), where=colsum > 0)
    return response, colsum


def dagostini_unfold(data, migration, prior, n_iter=10, data_var=None):
    """Unfold a reco distribution to true slots.

    Parameters
    ----------
    data      : (n_reco,) measured (background-subtracted) reco distribution.
    migration : (n_reco, n_true) M[reco, true], count-conserving.
    prior     : (n_true,) initial true-distribution guess (e.g. MC truth).
    n_iter    : D'Agostini iterations (paper: 10).
    data_var  : (n_reco,) optional per-bin variance of `data` for approximate
                error propagation.

    Returns
    -------
    unfolded  : (n_true,) unfolded true distribution.
    var       : (n_true,) propagated variance, or None if data_var is None.
    U         : (n_true, n_reco) final Bayes unfolding matrix.
    """
    response, colsum = response_matrix(migration)
    data = np.asarray(data, dtype=np.float64)
    n = np.asarray(prior, dtype=np.float64).copy()
    n[colsum == 0] = 0.0                                 # no info on empty true slots
    n_true, n_reco = migration.shape[1], migration.shape[0]
    U = np.zeros((n_true, n_reco))
    for _ in range(n_iter):
        folded = response @ n                            # (n_reco,) forward fold
        with np.errstate(invalid="ignore", divide="ignore"):
            U = np.where(folded > 0, (response * n).T / folded, 0.0)  # P(true|reco)
        n = U @ data
    var = None if data_var is None else (U ** 2) @ np.asarray(data_var, dtype=np.float64)
    return n, var, U
