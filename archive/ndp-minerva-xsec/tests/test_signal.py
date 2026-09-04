"""Tests for xsec.signal: boundary semantics + golden streamed counts.

Boundary tests are synthetic and offline. The golden test streams the Truth
tree of MC run 110040 via xrootd (marker: network) and must reproduce the
frozen-manifest anchors exactly: 397,604 truth-signal rows, of which 65,712
in phase space.
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from xsec import signal
from xsec.constants import (MAX_MU_THETA_RAD, NUMI_BEAM_ANGLE_RAD,
                            PZ_MIN_TRUE_MEV, Z_MAX_MM, Z_MIN_MM)

REPO_ROOT = Path(__file__).resolve().parents[1]

GOLDEN_TRUTH_SIGNAL = 397604   # frozen manifest + streamed 2026-06-12
GOLDEN_EFF_DENOM = 65712       # frozen manifest (runs/2026-06-09_dsigma_dpt)


def det_from_beam(xb, yb, zb):
    """Rotate a beam-frame 3-vector into detector coordinates.

    Inverse of the forward convention zb = sin(a)*py + cos(a)*pz used in
    xsec.kinematics.true_theta_p.
    """
    c, s = math.cos(NUMI_BEAM_ANGLE_RAD), math.sin(NUMI_BEAM_ANGLE_RAD)
    return xb, c * yb + s * zb, -s * yb + c * zb


def lepton(theta_rad, p_mev):
    """Detector-frame (px,py,pz,E) for a lepton at beam-angle theta."""
    xb, yb, zb = p_mev * math.sin(theta_rad), 0.0, p_mev * math.cos(theta_rad)
    px, py, pz = det_from_beam(xb, yb, zb)
    return [px, py, pz, p_mev]  # E column unused


GOOD_VTX = [0.0, 0.0, 7000.0, 0.0]
GOOD_LEP = lepton(0.1, 5000.0)


def ps(vtx_rows, lep_rows):
    return signal.in_phase_space(np.array(vtx_rows, dtype=float),
                                 np.array(lep_rows, dtype=float))


# --------------------------------------------------------------- is_signal
def test_is_signal_truth_table():
    inc = np.array([14, 14, -14, 12, 14])
    cur = np.array([1, 2, 1, 1, 0])
    assert signal.is_signal(inc, cur).tolist() == [True, False, False, False, False]


def test_signal_branch_lists_from_catalog():
    assert signal.SIGNAL_BRANCHES == ("mc_incoming", "mc_current")
    assert signal.PHASE_SPACE_BRANCHES == ("mc_vtx", "mc_primFSLepton")


# ----------------------------------------------------------- phase space: z
def test_z_boundaries_inclusive():
    vtxs = [[0, 0, Z_MIN_MM, 0], [0, 0, Z_MIN_MM - 1e-3, 0],
            [0, 0, Z_MAX_MM, 0], [0, 0, Z_MAX_MM + 1e-3, 0]]
    leps = [GOOD_LEP] * 4
    assert ps(vtxs, leps).tolist() == [True, False, True, False]


# --------------------------------------------------------- phase space: hex
def test_apothem_strict_boundaries():
    vtxs = [[0, 0, 7000, 0],
            [850.0, 0, 7000, 0],            # |x| == apothem -> out (strict)
            [849.999, 0, 7000, 0],
            [0, 981.0, 7000, 0],            # below 2*850/sqrt(3) ~ 981.495
            [0, 981.5, 7000, 0],            # above -> out
            [500.0, 692.0, 7000, 0],        # diagonal: limit ~ 692.82
            [500.0, 693.0, 7000, 0]]
    leps = [GOOD_LEP] * 7
    assert ps(vtxs, leps).tolist() == [True, False, True, True, False, True, False]


# ------------------------------------------------------- phase space: theta
def test_theta_boundary_inclusive():
    eps = 1e-6
    leps = [lepton(MAX_MU_THETA_RAD - eps, 5000.0),
            lepton(MAX_MU_THETA_RAD + eps, 5000.0)]
    vtxs = [GOOD_VTX] * 2
    assert ps(vtxs, leps).tolist() == [True, False]


# ---------------------------------------------------------- phase space: pz
def test_pz_boundary_inclusive():
    theta = 0.1
    p_at_limit = PZ_MIN_TRUE_MEV / math.cos(theta)
    leps = [lepton(theta, p_at_limit * (1 + 1e-9)),
            lepton(theta, p_at_limit * (1 - 1e-6))]
    vtxs = [GOOD_VTX] * 2
    assert ps(vtxs, leps).tolist() == [True, False]


def test_zero_momentum_fails():
    assert ps([GOOD_VTX], [[0.0, 0.0, 0.0, 0.0]]).tolist() == [False]


def test_efficiency_denominator_combines_both():
    inc, cur = np.array([14, 14]), np.array([1, 1])
    vtx = np.array([GOOD_VTX, [0, 0, 100.0, 0]])   # second outside fiducial
    lep = np.array([GOOD_LEP, GOOD_LEP])
    assert signal.is_efficiency_denominator(inc, cur, vtx, lep).tolist() == [True, False]


# ------------------------------------------------------------- golden gate
@pytest.mark.network
def test_golden_truth_counts_streamed():
    uproot = pytest.importorskip("uproot")
    spec = json.loads((REPO_ROOT / "config" / "datasets" / "me1A_single_pair.json").read_text())
    url = next(e["url"] for e in spec["files"] if e["role"] == "mc")
    try:
        f = uproot.open(url)
    except Exception as err:  # network/door unavailable -> skip, not fail
        pytest.skip(f"xrootd unreachable: {err}")
    with f:
        arrs = f["Truth"].arrays(["mc_incoming", "mc_current", "mc_vtx",
                                  "mc_primFSLepton"], library="np")
    sig = signal.is_signal(arrs["mc_incoming"], arrs["mc_current"])
    assert int(sig.sum()) == GOLDEN_TRUTH_SIGNAL
    denom = sig & signal.in_phase_space(arrs["mc_vtx"], arrs["mc_primFSLepton"])
    assert int(denom.sum()) == GOLDEN_EFF_DENOM
