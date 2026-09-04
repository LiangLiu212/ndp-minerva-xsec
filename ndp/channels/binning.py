"""Linearised 2D binning shared by channels, surrogates and comparisons.

A `Binning` has two axes (x = first observable, y = second) and a global-cell formula.
The two conventions in the wild are named exactly as the MINERvA benchmark manifests
name them so a channel can adopt a paper's linearisation verbatim:

    "ix*n_y + iy"   (2106.16210: GlobalID = ipt*n_pz + ipz, p_par inner)
    "iy*n_x + ix"   (2002.12496: GlobalID = ipz*n_pt + ipt, pT inner)

`axes` names the observables; edges are in the observable's units (GeV).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

FORMULAS = ("ix*n_y + iy", "iy*n_x + ix")
# Aliases so channel files can use the paper's variable names.
_ALIASES = {"ipt*n_pz + ipz": "ix*n_y + iy", "ipz*n_pt + ipt": "iy*n_x + ix"}


@dataclass(frozen=True)
class Binning:
    x_name: str
    y_name: str
    x_edges: tuple
    y_edges: tuple
    formula: str = "ix*n_y + iy"

    def __post_init__(self):
        object.__setattr__(self, "formula", _ALIASES.get(self.formula, self.formula))
        if self.formula not in FORMULAS:
            raise ValueError(f"unknown global-cell formula {self.formula!r}; known {FORMULAS} (+aliases {list(_ALIASES)})")
        for e in (self.x_edges, self.y_edges):
            if len(e) < 2 or np.any(np.diff(e) <= 0):
                raise ValueError("edges must be strictly increasing with >= 2 entries")

    @property
    def n_x(self) -> int:
        return len(self.x_edges) - 1

    @property
    def n_y(self) -> int:
        return len(self.y_edges) - 1

    @property
    def n_cells(self) -> int:
        return self.n_x * self.n_y

    def cell(self, ix: np.ndarray, iy: np.ndarray) -> np.ndarray:
        ix, iy = np.asarray(ix), np.asarray(iy)
        if self.formula == "ix*n_y + iy":
            return ix * self.n_y + iy
        return iy * self.n_x + ix

    def unravel(self, g: np.ndarray | None = None):
        """cell -> (ix, iy)."""
        g = np.arange(self.n_cells) if g is None else np.asarray(g)
        if self.formula == "ix*n_y + iy":
            return g // self.n_y, g % self.n_y
        return g % self.n_x, g // self.n_x

    def areas(self) -> np.ndarray:
        ix, iy = self.unravel()
        return np.diff(self.x_edges)[ix] * np.diff(self.y_edges)[iy]

    def digitize(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Global cell per event, or -1 when outside the grid (either axis)."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        xe, ye = np.asarray(self.x_edges), np.asarray(self.y_edges)
        ix = np.searchsorted(xe, x, side="right") - 1
        iy = np.searchsorted(ye, y, side="right") - 1
        inside = (x >= xe[0]) & (x < xe[-1]) & (y >= ye[0]) & (y < ye[-1])
        g = np.full(x.shape, -1, dtype=np.int64)
        g[inside] = self.cell(ix[inside], iy[inside])
        return g

    def histogram(self, x, y, weights=None):
        """(sumw[n_cells], sumw2[n_cells], n_outside) for events with observables x, y."""
        g = self.digitize(x, y)
        inside = g >= 0
        w = np.ones(len(g)) if weights is None else np.asarray(weights, float)
        sumw = np.bincount(g[inside], weights=w[inside], minlength=self.n_cells)
        sumw2 = np.bincount(g[inside], weights=w[inside] ** 2, minlength=self.n_cells)
        return sumw, sumw2, int((~inside).sum())

    def to_grid(self, cells: np.ndarray) -> np.ndarray:
        """cells[n_cells] -> 2D array [n_x, n_y]."""
        grid = np.zeros((self.n_x, self.n_y))
        ix, iy = self.unravel()
        grid[ix, iy] = np.asarray(cells)
        return grid

    def from_grid(self, grid: np.ndarray) -> np.ndarray:
        ix, iy = self.unravel()
        return np.asarray(grid)[ix, iy]

    def project(self, cells: np.ndarray, axis: str, per_width: bool = True) -> np.ndarray:
        """Sum cells over the other axis; divide by the projected bin width if per_width."""
        grid = self.to_grid(cells)
        if axis == "x":
            proj, w = grid.sum(axis=1), np.diff(self.x_edges)
        elif axis == "y":
            proj, w = grid.sum(axis=0), np.diff(self.y_edges)
        else:
            raise ValueError("axis must be 'x' or 'y'")
        return proj / w if per_width else proj

    def to_dict(self) -> dict:
        return {"x_name": self.x_name, "y_name": self.y_name, "x_edges": list(self.x_edges),
                "y_edges": list(self.y_edges), "formula": self.formula}

    @staticmethod
    def from_dict(d: dict) -> "Binning":
        return Binning(d["x_name"], d["y_name"], tuple(float(v) for v in d["x_edges"]),
                       tuple(float(v) for v in d["y_edges"]), d.get("formula", "ix*n_y + iy"))
