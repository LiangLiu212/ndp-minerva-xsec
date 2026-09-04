"""Parametric smearing surrogate: reco = f(true) + noise, learned per true cell, plus acceptance.

Where the binned response is a lookup table on a fixed grid, this model represents the
detector as a *conditional distribution* of the reconstructed observables given the true
ones. Per observable and per (coarse) true cell it stores the mean and standard deviation
of the reco/true ratio (or difference), and per true cell the acceptance (selection
efficiency). Sampling from it turns any truth-level event sample — on any binning — into
a reconstructed one, which is what "smearing generator truth to detector reconstruction"
means operationally. It is deliberately the simplest member of the family; the interface
(`fit` from paired MC, `sample_reco`, `fold_events`) is what a conditional normalising
flow or diffusion surrogate must implement to slot in.
"""
from __future__ import annotations

import numpy as np

from .base import Surrogate
from ..channels.binning import Binning


class SmearingSurrogate(Surrogate):
    kind = "parametric_smearing"

    def __init__(self, binning: Binning, eff: np.ndarray, ratio_mean: np.ndarray, ratio_std: np.ndarray,
                 den_counts: np.ndarray | None = None, mode: str = "diff,ratio", n_samples: int = 20,
                 meta: dict | None = None):
        super().__init__(binning, meta)
        n = binning.n_cells
        self.eff = np.asarray(eff, float).reshape(n)
        self.ratio_mean = np.asarray(ratio_mean, float).reshape(n, 2)   # per true cell, per observable (x, y)
        self.ratio_std = np.asarray(ratio_std, float).reshape(n, 2)
        self.den_counts = None if den_counts is None else np.asarray(den_counts, float).reshape(n)
        self.mode = mode                      # "<x mode>,<y mode>", each diff (reco = true + d) or ratio (reco = true * r)
        self.modes = tuple(mode.split(",")) if "," in mode else (mode, mode)
        self.n_samples = int(n_samples)

    # ---- the surrogate ---------------------------------------------------------------------
    def sample_reco(self, x_true, y_true, rng: np.random.Generator):
        """One reconstructed (x, y) per input event, or NaN when the event is not accepted."""
        g = self.binning.digitize(x_true, y_true)
        inside = g >= 0
        gi = np.where(inside, g, 0)
        accept = inside & (rng.random(len(g)) < self.eff[gi])
        mu, sd = self.ratio_mean[gi], self.ratio_std[gi]
        noise = rng.standard_normal((len(g), 2))
        out = []
        for k, tv in enumerate((x_true, y_true)):
            s_k = mu[:, k] + sd[:, k] * noise[:, k]
            out.append(tv * s_k if self.modes[k] == "ratio" else tv + s_k)
        xr, yr = out
        xr = np.where(accept, xr, np.nan)
        yr = np.where(accept, yr, np.nan)
        return xr, yr, accept

    def fold_events(self, x_true, y_true, weights=None, rng=None) -> np.ndarray:
        rng = np.random.default_rng(0) if rng is None else rng
        w = np.ones(len(x_true)) if weights is None else np.asarray(weights, float)
        acc = np.zeros(self.binning.n_cells)
        for _ in range(self.n_samples):
            xr, yr, ok = self.sample_reco(np.asarray(x_true, float), np.asarray(y_true, float), rng)
            h, _, _ = self.binning.histogram(xr[ok], yr[ok], w[ok])
            acc += h
        return acc / self.n_samples

    def fold(self, true_cells: np.ndarray) -> np.ndarray:
        """Fold an expected true-cell population by placing events at cell centres (log-uniform
        within the cell along each axis) — exact only in the limit of fine cells; prefer
        fold_events when the truth events are available."""
        rng = np.random.default_rng(1)
        true_cells = np.asarray(true_cells, float)
        ix, iy = self.binning.unravel()
        xe, ye = np.asarray(self.binning.x_edges), np.asarray(self.binning.y_edges)
        reps = 200
        xs = np.repeat(xe[ix], reps) + np.tile(rng.random(reps), self.binning.n_cells) * np.repeat(np.diff(xe)[ix], reps)
        ys = np.repeat(ye[iy], reps) + np.tile(rng.random(reps), self.binning.n_cells) * np.repeat(np.diff(ye)[iy], reps)
        ws = np.repeat(true_cells / reps, reps)
        return self.fold_events(xs, ys, ws, rng)

    def diagnostics(self) -> dict:
        have = self.den_counts is not None and (self.den_counts > 0).any()
        return {"n_cells": self.binning.n_cells, "modes": list(self.modes), "n_samples": self.n_samples,
                "eff_population_weighted_mean": float(np.sum(self.eff * self.den_counts) / self.den_counts.sum()) if have else None,
                "median_width_x": float(np.nanmedian(self.ratio_std[:, 0])),
                "median_width_y": float(np.nanmedian(self.ratio_std[:, 1])),
                "median_location_x": float(np.nanmedian(self.ratio_mean[:, 0])),
                "median_location_y": float(np.nanmedian(self.ratio_mean[:, 1]))}

    # ---- construction ------------------------------------------------------------------------
    @classmethod
    def fit(cls, binning: Binning, *, x_true_den, y_true_den, x_true_num, y_true_num, x_reco_num, y_reco_num,
            w_den=None, w_num=None, mode: str = "diff,ratio", min_events: int = 20, n_samples: int = 20,
            robust: bool = True, meta: dict | None = None) -> "SmearingSurrogate":
        """Per true cell: location + width of (reco - true) [diff] or reco/true [ratio] per axis.

        `robust=True` uses median and 1.4826*MAD (the Gaussian-core resolution, insensitive
        to the non-Gaussian tails a mean/std would be dominated by); the tails are then NOT
        modelled — that is the documented limitation of this surrogate class."""
        n = binning.n_cells
        den, _, _ = binning.histogram(x_true_den, y_true_den, w_den)
        num, _, _ = binning.histogram(x_true_num, y_true_num, w_num)
        with np.errstate(invalid="ignore", divide="ignore"):
            eff = np.clip(np.where(den > 0, num / den, 0.0), 0.0, 1.0)
        gt = binning.digitize(x_true_num, y_true_num)
        ok = gt >= 0
        xt, yt, xr, yr = (np.asarray(a, float)[ok] for a in (x_true_num, y_true_num, x_reco_num, y_reco_num))
        g = gt[ok]
        modes = tuple(mode.split(",")) if "," in mode else (mode, mode)
        res = []
        for k, (tv, rv) in enumerate(((xt, xr), (yt, yr))):
            res.append(rv / np.where(tv != 0, tv, np.nan) if modes[k] == "ratio" else rv - tv)
        rx, ry = res
        mean = np.zeros((n, 2)); std = np.zeros((n, 2)); cnt = np.bincount(g, minlength=n)
        order = np.argsort(g, kind="stable")
        g_sorted = g[order]
        bounds = np.searchsorted(g_sorted, np.arange(n + 1))
        for k, r in enumerate((rx, ry)):
            r_sorted = r[order]
            for j in range(n):
                seg = r_sorted[bounds[j]:bounds[j + 1]]
                seg = seg[np.isfinite(seg)]
                if seg.size == 0:
                    mean[j, k], std[j, k] = np.nan, np.nan
                elif robust:
                    med = np.median(seg)
                    mean[j, k], std[j, k] = med, 1.4826 * np.median(np.abs(seg - med))
                else:
                    mean[j, k], std[j, k] = seg.mean(), seg.std(ddof=1) if seg.size > 1 else np.nan
        # cells with too few events borrow the global resolution
        fx, fy = np.isfinite(rx), np.isfinite(ry)
        if robust:
            glob_mean = np.array([np.median(rx[fx]), np.median(ry[fy])])
            glob_std = np.array([1.4826 * np.median(np.abs(rx[fx] - glob_mean[0])), 1.4826 * np.median(np.abs(ry[fy] - glob_mean[1]))])
        else:
            glob_mean = np.array([np.mean(rx[fx]), np.mean(ry[fy])]); glob_std = np.array([np.std(rx[fx]), np.std(ry[fy])])
        sparse = cnt < min_events
        mean[sparse] = glob_mean; std[sparse] = glob_std
        mean = np.where(np.isfinite(mean), mean, glob_mean); std = np.where(np.isfinite(std), std, glob_std)
        m = dict(meta or {})
        m.update({"n_den": float(den.sum()), "n_num": float(num.sum()), "n_sparse_cells": int(sparse.sum()),
                  "min_events": min_events, "robust": robust, "modes": list(modes),
                  "global_location": glob_mean.tolist(), "global_width": glob_std.tolist()})
        return cls(binning, eff, mean, std, den_counts=den, mode=mode, n_samples=n_samples, meta=m)

    def _arrays(self) -> dict:
        d = {"eff": self.eff, "ratio_mean": self.ratio_mean, "ratio_std": self.ratio_std,
             "mode": np.array(self.mode), "n_samples": np.array(self.n_samples)}
        if self.den_counts is not None:
            d["den_counts"] = self.den_counts
        return d

    @classmethod
    def _from_arrays(cls, binning, arrays, meta):
        return cls(binning, arrays["eff"], arrays["ratio_mean"], arrays["ratio_std"], den_counts=arrays.get("den_counts"),
                   mode=str(arrays["mode"]), n_samples=int(arrays["n_samples"]), meta=meta)
