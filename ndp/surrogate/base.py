"""Detector-surrogate interface: p(reco | truth) for one channel's observables.

A surrogate answers one question — given a true-cell (or true-event) population of the
channel's signal, what reconstructed, selected population does the detector produce —
and carries its own provenance (what MC it was learned from, when, how much).

    fold(true_cells)          expected reco cells from an expected true-cell population
    fold_events(table, w)     same, from truth events + weights (uses per-event smearing
                              when the surrogate has it, else bins first)
    background(pot_data)      expected non-signal reco counts at the given exposure (or 0)
"""
from __future__ import annotations

import json
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np

from ..channels.binning import Binning


class Surrogate(ABC):
    kind: str = "abstract"

    def __init__(self, binning: Binning, meta: dict | None = None):
        self.binning = binning
        self.meta = dict(meta or {})

    @abstractmethod
    def fold(self, true_cells: np.ndarray) -> np.ndarray: ...

    def fold_events(self, x_true, y_true, weights, rng=None) -> np.ndarray:
        cells, _, _ = self.binning.histogram(x_true, y_true, weights)
        return self.fold(cells)

    def background(self, pot_data: float) -> np.ndarray:
        return np.zeros(self.binning.n_cells)

    def fold_variance(self, true_cells: np.ndarray) -> np.ndarray:
        """MC-statistics variance of fold(true_cells); zero unless a subclass estimates it."""
        return np.zeros(self.binning.n_cells)

    # ---- persistence -----------------------------------------------------------------
    def _arrays(self) -> dict:
        return {}

    def save(self, directory: str | Path) -> Path:
        d = Path(directory)
        d.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(d / "arrays.npz", **self._arrays())
        (d / "surrogate.json").write_text(json.dumps(
            {"kind": self.kind, "binning": self.binning.to_dict(), "meta": self.meta}, indent=2, default=str))
        return d

    @classmethod
    def _from_arrays(cls, binning: Binning, arrays: dict, meta: dict) -> "Surrogate":
        raise NotImplementedError


def load_surrogate(directory: str | Path) -> Surrogate:
    from . import binned, parametric
    d = Path(directory)
    spec = json.loads((d / "surrogate.json").read_text())
    z = np.load(d / "arrays.npz", allow_pickle=False)
    arrays = {k: z[k] for k in z.files}
    binning = Binning.from_dict(spec["binning"])
    kinds = {"binned_response": binned.BinnedResponse, "parametric_smearing": parametric.SmearingSurrogate}
    return kinds[spec["kind"]]._from_arrays(binning, arrays, spec.get("meta", {}))
