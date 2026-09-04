"""Binned detector response: efficiency(true cell) x migration P(reco cell | true cell) + background.

This is the classical forward-folding surrogate. It is exact for the binning it was
built on and needs no assumption about the smearing shape; its price is that it cannot
be evaluated on a different binning and it inherits the MC statistics per cell.

    reco[i] = sum_j P[i, j] * eff[j] * true[j] + bkg[i]

Built either from a MINERvA MC AnaTuple (reco tree for numerator + migration, Truth tree
for the denominator) or from the exported artifacts of an audited run.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .base import Surrogate
from ..channels.binning import Binning


class BinnedResponse(Surrogate):
    kind = "binned_response"

    def __init__(self, binning: Binning, eff: np.ndarray, P: np.ndarray, migration_counts: np.ndarray | None = None,
                 den_counts: np.ndarray | None = None, num_counts: np.ndarray | None = None,
                 bkg_per_pot: np.ndarray | None = None, meta: dict | None = None):
        super().__init__(binning, meta)
        n = binning.n_cells
        self.eff = np.asarray(eff, float).reshape(n)
        self.P = np.asarray(P, float).reshape(n, n)
        self.migration_counts = None if migration_counts is None else np.asarray(migration_counts, float).reshape(n, n)
        self.den_counts = None if den_counts is None else np.asarray(den_counts, float).reshape(n)
        self.num_counts = None if num_counts is None else np.asarray(num_counts, float).reshape(n)
        self.bkg_per_pot = None if bkg_per_pot is None else np.asarray(bkg_per_pot, float).reshape(n)

    # ---- the surrogate -------------------------------------------------------------------
    def fold(self, true_cells: np.ndarray) -> np.ndarray:
        true_cells = np.asarray(true_cells, float).reshape(self.binning.n_cells)
        return self.P @ (self.eff * true_cells)

    def background(self, pot_data: float) -> np.ndarray:
        if self.bkg_per_pot is None:
            return np.zeros(self.binning.n_cells)
        return self.bkg_per_pot * pot_data

    def fold_variance(self, true_cells: np.ndarray) -> np.ndarray:
        """Approximate MC-stat variance: binomial on eff (per true cell) + multinomial on P
        columns, propagated linearly and treated as independent across true cells."""
        if self.den_counts is None or self.migration_counts is None:
            return np.zeros(self.binning.n_cells)
        t = np.asarray(true_cells, float)
        with np.errstate(invalid="ignore", divide="ignore"):
            var_eff = np.where(self.den_counts > 0, self.eff * (1 - self.eff) / self.den_counts, 0.0)
            col = self.migration_counts.sum(axis=0)
            var_P = np.where(col > 0, self.P * (1 - self.P) / col, 0.0)
        var = (self.P ** 2) @ (var_eff * t ** 2) + var_P @ ((self.eff * t) ** 2)
        return var

    def eff_uncertainty(self) -> np.ndarray:
        if self.den_counts is None:
            return np.zeros_like(self.eff)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(self.den_counts > 0, np.sqrt(self.eff * (1 - self.eff) / self.den_counts), 0.0)

    def diagnostics(self) -> dict:
        col = self.P.sum(axis=0)
        diag = np.diag(self.P)
        out = {"n_cells": self.binning.n_cells,
               "eff_min": float(np.nanmin(self.eff)), "eff_max": float(np.nanmax(self.eff)),
               "n_true_cells_with_response": int((col > 0).sum()),
               "P_column_sums_range": [float(col[col > 0].min()) if (col > 0).any() else None, float(col.max())],
               "median_diagonal_P": float(np.median(diag[col > 0])) if (col > 0).any() else None}
        if self.den_counts is not None:
            w = self.den_counts
            out["eff_population_weighted_mean"] = float(np.sum(self.eff * w) / w.sum()) if w.sum() else None
            out["n_den"] = float(w.sum())
        if self.num_counts is not None:
            out["n_num"] = float(self.num_counts.sum())
        if self.migration_counts is not None:
            out["n_migration"] = float(self.migration_counts.sum())
        return out

    # ---- construction ----------------------------------------------------------------------
    @classmethod
    def fit(cls, binning: Binning, *, x_true_den, y_true_den, w_den=None,
            x_true_num, y_true_num, x_reco_num, y_reco_num, w_num=None,
            x_reco_bkg=None, y_reco_bkg=None, w_bkg=None, pot_mc: float | None = None,
            meta: dict | None = None) -> "BinnedResponse":
        """Learn eff, P and the additive background from paired MC.

        den : all signal events in the true phase space (Truth tree) -> true cells.
        num : reco-selected signal events in the true phase space, with true AND reco
              observables. eff[j] = num[j]/den[j]. P[:, j] = M[:, j] / num[j] where M
              counts pairs with both cells inside the grid — so a column sums to < 1 when
              part of that true cell reconstructs outside the reco grid (a real loss,
              not redistributed).
        bkg : every other reco-selected event inside the reco grid — non-signal AND
              signal whose truth is outside the phase space / true grid (the feed-in that
              a phase-space-restricted response cannot produce). Stored per POT.

        With these definitions folding the training MC's own truth reproduces its
        reco-selected in-grid count exactly (closure by construction).
        """
        n = binning.n_cells
        den, _, _ = binning.histogram(x_true_den, y_true_den, w_den)
        gt = binning.digitize(x_true_num, y_true_num)
        gr = binning.digitize(x_reco_num, y_reco_num)
        w = np.ones(len(gt)) if w_num is None else np.asarray(w_num, float)
        in_true = gt >= 0
        num = np.bincount(gt[in_true], weights=w[in_true], minlength=n)
        with np.errstate(invalid="ignore", divide="ignore"):
            eff = np.where(den > 0, num / den, 0.0)
        eff = np.clip(eff, 0.0, 1.0)
        both = in_true & (gr >= 0)
        M = np.zeros((n, n))
        np.add.at(M, (gr[both], gt[both]), w[both])
        P = np.zeros_like(M)
        P[:, num > 0] = M[:, num > 0] / num[num > 0]
        # feed-in: numerator events whose TRUE cell is outside the grid but reco inside
        feed = (~in_true) & (gr >= 0)
        bkg_counts = np.bincount(gr[feed], weights=w[feed], minlength=n)
        n_feed_true_out = float(w[feed].sum())
        if x_reco_bkg is not None:
            b, _, _ = binning.histogram(x_reco_bkg, y_reco_bkg, w_bkg)
            bkg_counts = bkg_counts + b
        bkg = bkg_counts / pot_mc if pot_mc else None
        m = dict(meta or {})
        m.update({"n_den": float(den.sum()), "n_num": float(num.sum()), "n_migration": float(M.sum()),
                  "n_num_reco_out_of_grid": float(w[in_true & (gr < 0)].sum()),
                  "n_feedin_true_out_of_grid": n_feed_true_out,
                  "n_bkg_total_in_reco_grid": float(bkg_counts.sum()), "pot_mc": pot_mc,
                  "response_convention": "P[reco, true] = M/num (column sums <= 1); bkg = non-signal + out-of-phase-space signal, per POT"})
        return cls(binning, eff, P, migration_counts=M, den_counts=den, num_counts=num, bkg_per_pot=bkg, meta=m)

    @classmethod
    def from_run_artifacts(cls, binning: Binning, migration_npy: str | Path, efficiency_npy: str | Path,
                           artifact_formula: str = "iy*n_x + ix", meta: dict | None = None) -> "BinnedResponse":
        """Build from an audited run's exported migration (reco x true) and efficiency grid.

        The MINERvA run 2026-06-19 stored migration in its *internal* linearisation
        g = ipl*14 + ipt (p_par outer, pT inner) and efficiency as a [16, 14] (p_par, pT)
        grid. Both are re-indexed into this binning's convention here.
        """
        M_int = np.load(migration_npy)
        eff_grid = np.load(efficiency_npy)          # [n_y(pz)=16, n_x(pt)=14]
        n = binning.n_cells
        if M_int.shape != (n, n):
            raise ValueError(f"migration shape {M_int.shape} != ({n},{n})")
        src = Binning(binning.x_name, binning.y_name, binning.x_edges, binning.y_edges, artifact_formula)
        # permutation: for every (ix, iy) map src cell -> dst cell
        ix, iy = src.unravel()
        perm = binning.cell(ix, iy)                  # dst index for each src index
        M = np.zeros_like(M_int)
        M[np.ix_(perm, perm)] = M_int
        eff_cells = np.zeros(n)
        eff_grid = np.nan_to_num(np.asarray(eff_grid, float))
        eff_cells[binning.cell(*np.meshgrid(np.arange(binning.n_x), np.arange(binning.n_y), indexing="ij"))] = eff_grid.T
        col = M.sum(axis=0)
        P = np.zeros_like(M)
        P[:, col > 0] = M[:, col > 0] / col[col > 0]
        m = dict(meta or {})
        m.update({"built_from": {"migration": str(migration_npy), "efficiency": str(efficiency_npy),
                                 "artifact_formula": artifact_formula},
                  "n_migration": float(M.sum()), "note": "denominator counts not exported by the run; "
                  "efficiency uncertainties unavailable"})
        return cls(binning, eff_cells, P, migration_counts=M, meta=m)

    # ---- persistence ---------------------------------------------------------------------
    def _arrays(self) -> dict:
        d = {"eff": self.eff, "P": self.P}
        for k in ("migration_counts", "den_counts", "num_counts", "bkg_per_pot"):
            v = getattr(self, k)
            if v is not None:
                d[k] = v
        return d

    @classmethod
    def _from_arrays(cls, binning, arrays, meta):
        return cls(binning, arrays["eff"], arrays["P"], migration_counts=arrays.get("migration_counts"),
                   den_counts=arrays.get("den_counts"), num_counts=arrays.get("num_counts"),
                   bkg_per_pot=arrays.get("bkg_per_pot"), meta=meta)
