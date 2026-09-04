"""Tests for xsec.cuts: per-cut boundaries, composition, golden parity.

Boundary tests are synthetic and offline. The golden tests (marker: network)
stream the selection branches of the golden pair and must reproduce the
certified counts exactly: 844 selected data / 43643 selected MC, splitting
into 43539 signal + 104 background (47 NC + 57 Wrong_Sign).
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from xsec import cuts, signal
from xsec.constants import MAX_MU_THETA_RAD, NC_CURRENT, Z_MAX_MM, Z_MIN_MM

REPO_ROOT = Path(__file__).resolve().parents[1]

GOLDEN_DATA_SELECTED = 844
GOLDEN_MC_SELECTED = 43643
GOLDEN_MC_SIGNAL = 43539
GOLDEN_MC_BACKGROUND = 104
GOLDEN_MC_BKG_NC = 47
GOLDEN_MC_BKG_WRONG_SIGN = 57


def good_event():
    """Branch dict for one event that passes all six cuts."""
    return {
        "vtx": np.array([[0.0, 0.0, 7000.0, 0.0]]),
        "muon_thetaX": np.array([0.05]),
        "muon_thetaY": np.array([0.05]),
        "isMinosMatchTrack": np.array([1], dtype=np.int32),
        "phys_n_dead_discr_pair_upstream_prim_track_proj": np.array([0], dtype=np.int32),
        "MasterAnaDev_minos_trk_qp": np.array([-0.001]),
    }


def with_(key, value):
    ev = good_event()
    ev[key] = np.asarray(value) if key != "vtx" else np.array(value)
    return ev


def test_good_event_passes():
    assert cuts.reco_selection(good_event()).tolist() == [True]


def test_branch_list_from_catalog():
    assert set(cuts.RECO_SELECTION_BRANCHES) == set(good_event().keys())


def test_z_range_boundaries_inclusive():
    for z, expect in [(Z_MIN_MM, True), (Z_MIN_MM - 1e-3, False),
                      (Z_MAX_MM, True), (Z_MAX_MM + 1e-3, False)]:
        ev = with_("vtx", [[0.0, 0.0, z, 0.0]])
        assert cuts.reco_selection(ev).tolist() == [expect], f"z={z}"


def test_apothem_strict():
    for xy, expect in [((850.0, 0.0), False), ((849.999, 0.0), True),
                       ((0.0, 981.0), True), ((0.0, 981.5), False),
                       ((500.0, 692.0), True), ((500.0, 693.0), False)]:
        ev = with_("vtx", [[xy[0], xy[1], 7000.0, 0.0]])
        assert cuts.reco_selection(ev).tolist() == [expect], f"xy={xy}"


def test_muon_angle_strict():
    eps = 1e-6
    for tx, expect in [(MAX_MU_THETA_RAD - eps, True),
                       (MAX_MU_THETA_RAD + eps, False)]:
        ev = good_event()
        ev["muon_thetaX"] = np.array([tx])
        ev["muon_thetaY"] = np.array([0.0])
        # theta3d(tx, 0) == tx
        assert cuts.reco_selection(ev).tolist() == [expect], f"thetaX={tx}"


def test_minos_match_exact_one():
    for v, expect in [(1, True), (0, False), (-1, False), (2, False)]:
        ev = with_("isMinosMatchTrack", np.array([v], dtype=np.int32))
        assert cuts.reco_selection(ev).tolist() == [expect], f"match={v}"


def test_deadtime_at_most_one():
    for v, expect in [(0, True), (1, True), (2, False)]:
        ev = with_("phys_n_dead_discr_pair_upstream_prim_track_proj",
                   np.array([v], dtype=np.int32))
        assert cuts.reco_selection(ev).tolist() == [expect], f"dead={v}"


def test_neutrino_negative_qp():
    for v, expect in [(-1e-6, True), (0.0, False), (1e-6, False)]:
        ev = with_("MasterAnaDev_minos_trk_qp", np.array([v]))
        assert cuts.reco_selection(ev).tolist() == [expect], f"qp={v}"


def test_unmatched_track_fails_regardless_of_qp():
    ev = with_("isMinosMatchTrack", np.array([0], dtype=np.int32))
    ev["MasterAnaDev_minos_trk_qp"] = np.array([-999.0])  # sentinel-ish
    assert cuts.reco_selection(ev).tolist() == [False]


def test_cutflow_order_and_final_count():
    evs = {k: np.concatenate([good_event()[k], with_("vtx", [[0, 0, 100.0, 0]])[k]])
           for k in good_event()}
    flow = cuts.cutflow(evs)
    assert [lbl for lbl, _ in flow] == list(cuts.CUT_LABELS)
    assert flow[0] == ("ZRange", 1)            # bad-z event dies first
    assert flow[-1][1] == int(cuts.reco_selection(evs).sum())


# ------------------------------------------------------------- golden gates
def _stream(role, branches):
    uproot = pytest.importorskip("uproot")
    spec = json.loads((REPO_ROOT / "config" / "datasets" / "me1A_single_pair.json").read_text())
    url = next(e["url"] for e in spec["files"] if e["role"] == role)
    try:
        f = uproot.open(url)
    except Exception as err:
        pytest.skip(f"xrootd unreachable: {err}")
    with f:
        return f["MasterAnaDev"].arrays(branches, library="np")


@pytest.mark.network
def test_golden_data_selection():
    arrs = _stream("data", list(cuts.RECO_SELECTION_BRANCHES))
    assert int(cuts.reco_selection(arrs).sum()) == GOLDEN_DATA_SELECTED


@pytest.mark.network
def test_golden_mc_selection_and_split():
    arrs = _stream("mc", list(cuts.RECO_SELECTION_BRANCHES)
                   + list(signal.SIGNAL_BRANCHES))
    selected = cuts.reco_selection(arrs)
    assert int(selected.sum()) == GOLDEN_MC_SELECTED

    sig = signal.is_signal(arrs["mc_incoming"], arrs["mc_current"])
    n_signal = int((selected & sig).sum())
    n_bkg = int((selected & ~sig).sum())
    assert n_signal == GOLDEN_MC_SIGNAL
    assert n_bkg == GOLDEN_MC_BACKGROUND

    nc = selected & ~sig & (np.asarray(arrs["mc_current"]) == NC_CURRENT)
    assert int(nc.sum()) == GOLDEN_MC_BKG_NC
    assert n_bkg - int(nc.sum()) == GOLDEN_MC_BKG_WRONG_SIGN
