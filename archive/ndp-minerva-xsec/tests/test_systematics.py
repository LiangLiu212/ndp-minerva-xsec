"""Tests for xsec.systematics: covariance algebra (CalcCovMx forms)."""
import numpy as np

from xsec import systematics as sx


def test_pair_covariance_outer_half_difference():
    plus = np.array([3.0, 1.0])
    minus = np.array([1.0, -1.0])
    cov = sx.pair_covariance(plus, minus)          # h = (1, 1)
    assert np.allclose(cov, [[1.0, 1.0], [1.0, 1.0]])


def test_sample_covariance_analytic():
    # 4 universes about mean 0: ±1 in bin0, ±2 in bin1
    U = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 2.0], [0.0, -2.0]])
    cov = sx.sample_covariance(U)
    assert np.allclose(cov, [[0.5, 0.0], [0.0, 2.0]])


def test_sample_covariance_matches_numpy_population_cov():
    rng = np.random.default_rng(3)
    U = rng.normal(size=(200, 5))
    cov = sx.sample_covariance(U)
    ref = np.cov(U, rowvar=False, bias=True)        # divide by N (not N-1)
    assert np.allclose(cov, ref, atol=1e-12)


def test_covariance_symmetric_and_psd():
    rng = np.random.default_rng(4)
    U = rng.normal(size=(50, 8))
    cov = sx.sample_covariance(U)
    assert np.allclose(cov, cov.T)
    assert np.linalg.eigvalsh(cov).min() > -1e-10   # PSD


def test_total_covariance_is_sum():
    a = np.array([[1.0, 0.0], [0.0, 1.0]])
    b = np.array([[2.0, 0.5], [0.5, 2.0]])
    stat = np.diag([0.1, 0.2])
    assert np.allclose(sx.total_covariance([a, b], stat), a + b + stat)
    assert np.allclose(sx.total_covariance([a, b]), a + b)


def test_normalization_covariance():
    sigma = np.array([10.0, 20.0])
    cov = sx.normalization_covariance(sigma, 0.04)      # 4% normalization
    # diagonal fractional error == frac; fully correlated
    assert np.allclose(sx.fractional_error(cov, sigma), 0.04)
    assert np.allclose(sx.correlation_matrix(cov), 1.0)


def test_fractional_error_and_correlation():
    cov = np.array([[4.0, 1.0], [1.0, 9.0]])
    frac = sx.fractional_error(cov, np.array([2.0, 3.0]))
    assert np.allclose(frac, [1.0, 1.0])            # √4/2, √9/3
    corr = sx.correlation_matrix(cov)
    assert np.allclose(np.diag(corr), 1.0)
    assert np.isclose(corr[0, 1], 1.0 / (2 * 3))
