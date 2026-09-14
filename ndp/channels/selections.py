"""Reconstruction-level selections, dispatched by a channel manifest's `selection.name`.

A selection is a function of the cached reco table (`ndp data cache`) returning a boolean mask
of selected candidates; `cutflow` returns the ordered cumulative steps. Cut values are read from
the manifest's `selection.params` — nothing is defaulted in code.

    minerva_cc_inclusive_v1     the MINERvA-101 chain, evaluated at cache-build time (`passed`,
                                `failing_cut`, certified row-by-row against the exploration repo's tool)
    minerva_ccqelike_1mu1p_v0   arXiv:2503.15047 transcription on the v3 cache columns
                                (adapters/minerva_anatuple.py::ccqelike_1mu1p_cutflow)
"""
from __future__ import annotations

import numpy as np

from ..adapters.minerva_anatuple import (CUT_LABELS, CCQELIKE_1MU1P_LABELS, ccqelike_1mu1p_cutflow,
                                         RECO_CACHE_VERSION)


def _n(r: dict) -> int:
    return len(next(iter(v for k, v in r.items() if k != "__meta__")))


def _cumulative(failing: np.ndarray, labels) -> list[tuple[str, np.ndarray]]:
    out = []
    for i, lab in enumerate(labels):
        out.append((lab, (failing == -1) | (failing > i)))
    return out


def _cc_inclusive(channel, r: dict):
    passed = np.asarray(r["passed"], bool)
    if "failing_cut" in r:
        return passed, np.asarray(r["failing_cut"], np.int64), CUT_LABELS
    return passed, np.where(passed, -1, 0), ("passed",)          # a table carrying only the final flag (toy / legacy)


def _ccqelike_1mu1p(channel, r: dict):
    meta = r.get("__meta__", {}) or {}
    if int(meta.get("cache_version", 1)) < 3 or "reco_proton_p" not in r:
        raise ValueError("selection minerva_ccqelike_1mu1p_v0 needs a version-3 reco cache "
                         f"(have {meta.get('cache_version', '?')}); rebuild with `python -m ndp data cache --channel {channel.name}`")
    params = channel.selection.get("params")
    if not params:
        raise KeyError(f"channel {channel.name}: selection.params missing (cut values live in the manifest)")
    passed, failing = ccqelike_1mu1p_cutflow(r, params)
    return passed, failing, CCQELIKE_1MU1P_LABELS


SELECTIONS = {"minerva_cc_inclusive_v1": _cc_inclusive, "minerva_ccqelike_1mu1p_v0": _ccqelike_1mu1p}


def _dispatch(channel):
    name = str(channel.selection.get("name", "minerva_cc_inclusive_v1"))
    if name not in SELECTIONS:
        raise KeyError(f"unknown selection {name!r}; known: {sorted(SELECTIONS)}")
    return SELECTIONS[name]


def select(channel, r: dict) -> np.ndarray:
    """Boolean mask of the reco candidates the channel's selection keeps."""
    return _dispatch(channel)(channel, r)[0]


def cutflow(channel, r: dict) -> list[tuple[str, np.ndarray]]:
    """Ordered (label, cumulative mask) steps of the channel's selection (all candidates first)."""
    passed, failing, labels = _dispatch(channel)(channel, r)
    return [("all_candidates", np.ones(_n(r), bool))] + _cumulative(failing, labels)
