"""Covariance assembly for the systematic-uncertainty stage (M4).

A systematic is an ensemble of universes, each a full cross-section re-run with
one input varied. The covariance is the sample covariance of the per-universe
cross sections, reproducing MnvVertErrorBand::CalcCovMx:

  many-universe (e.g. 100 flux PPFX):
      Cov[i,k] = (1/ΣW) Σ_u W_u (σ_u[i]-σ̄[i])(σ_u[k]-σ̄[k]),  σ̄ = weighted mean
  ±1σ pair (e.g. muon energy scale):
      Cov[i,k] = h[i] h[k],  h = (σ⁺ - σ⁻)/2

The total covariance is the sum of all group covariances plus the statistical
covariance (MnvH1D::GetTotalErrorMatrix). All functions are pure NumPy on
flat cross-section vectors (224 measurement cells).
"""
import numpy as np


def sample_covariance(universes, mean=None, weights=None):
    """Many-universe sample covariance (about the weighted universe mean).

    universes : (n_univ, n_bins) per-universe cross sections.
    mean      : optional (n_bins,) center; default = weighted universe mean.
    weights   : optional (n_univ,) universe weights; default uniform.
    """
    U = np.asarray(universes, dtype=np.float64)
    w = np.ones(U.shape[0]) if weights is None else np.asarray(weights, np.float64)
    W = w.sum()
    m = (w[:, None] * U).sum(0) / W if mean is None else np.asarray(mean, np.float64)
    d = U - m
    return np.einsum("u,ui,uj->ij", w, d, d) / W


def pair_covariance(plus, minus):
    """±1σ pair covariance = outer((σ⁺-σ⁻)/2) (CalcCovMx 2-universe form)."""
    h = (np.asarray(plus, np.float64) - np.asarray(minus, np.float64)) / 2.0
    return np.outer(h, h)


def total_covariance(group_covs, cov_stat=None):
    """Σ group covariances (+ optional statistical covariance)."""
    tot = np.zeros_like(np.asarray(next(iter(group_covs)), np.float64))
    for c in group_covs:
        tot = tot + np.asarray(c, np.float64)
    if cov_stat is not None:
        tot = tot + np.asarray(cov_stat, np.float64)
    return tot


def fractional_error(cov, sigma_cv):
    """Per-bin fractional uncertainty √diag(cov) / |σ_cv|."""
    d = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    sigma_cv = np.asarray(sigma_cv, np.float64)
    return np.divide(d, np.abs(sigma_cv), out=np.zeros_like(d), where=sigma_cv != 0)


def correlation_matrix(cov):
    """Correlation matrix from a covariance (0 where a diagonal is 0)."""
    d = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    denom = np.outer(d, d)
    return np.divide(cov, denom, out=np.zeros_like(cov, dtype=np.float64),
                     where=denom > 0)
