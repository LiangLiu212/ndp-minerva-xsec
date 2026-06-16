"""Tests for xsec.unfold: response normalization, closure, convergence."""
import numpy as np
import pytest

from xsec import unfold


def fold(migration, truth):
    """Reco distribution from folding a true vector through the migration."""
    R, _ = unfold.response_matrix(migration)
    return R @ truth


def test_response_columns_normalized():
    M = np.array([[8.0, 1.0], [2.0, 9.0]])     # [reco, true]
    R, colsum = unfold.response_matrix(M)
    assert np.allclose(R.sum(axis=0), 1.0)     # columns (true) sum to 1
    assert np.allclose(colsum, M.sum(axis=0))


def test_empty_true_column_stays_zero():
    M = np.array([[5.0, 0.0], [3.0, 0.0]])
    R, colsum = unfold.response_matrix(M)
    assert colsum[1] == 0 and np.all(R[:, 1] == 0)


def test_exact_closure_with_true_prior():
    # count-conserving migration; truth = column sums; data = folded truth.
    # Unfolding with prior == truth recovers truth to machine precision.
    M = np.array([[40.0, 5.0, 1.0],
                  [8.0, 50.0, 6.0],
                  [1.0, 7.0, 30.0]])
    truth = M.sum(axis=0)
    data = fold(M, truth)
    n, _, _ = unfold.dagostini_unfold(data, M, prior=truth, n_iter=1)
    assert np.allclose(n, truth, rtol=0, atol=1e-9)


def test_convergence_from_flat_prior():
    # exact folded data + count-conserving migration -> D'Agostini converges
    # to the truth regardless of prior; 50 iterations get there tightly.
    rng = np.random.default_rng(1)
    M = rng.uniform(0, 1, (5, 5)) + 5 * np.eye(5)   # diagonally dominant
    truth = M.sum(axis=0)
    data = fold(M, truth)
    flat = np.full(5, truth.mean())
    n, _, _ = unfold.dagostini_unfold(data, M, prior=flat, n_iter=50)
    assert np.allclose(n, truth, rtol=1e-6)


def test_two_bin_migration_recovers_shape():
    # asymmetric 2-bin migration, non-flat truth, exact prior
    M = np.array([[90.0, 10.0], [10.0, 90.0]])
    truth = np.array([70.0, 30.0])
    data = fold(M, truth)
    n, _, _ = unfold.dagostini_unfold(data, M, prior=truth, n_iter=4)
    assert np.allclose(n, truth, atol=1e-9)


def test_var_propagation_shape_and_positivity():
    M = np.array([[40.0, 5.0], [8.0, 50.0]])
    truth = M.sum(axis=0)
    data = fold(M, truth)
    n, var, U = unfold.dagostini_unfold(data, M, prior=truth, n_iter=10,
                                        data_var=data)
    assert var.shape == truth.shape
    assert np.all(var >= 0)
    assert U.shape == (M.shape[1], M.shape[0])
