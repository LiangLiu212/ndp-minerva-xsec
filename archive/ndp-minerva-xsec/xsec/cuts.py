"""The 6-cut CC-inclusive reco selection, vectorized.

Port of the certified ``tools/cc_inclusive_selector.py`` (exploration repo,
0 mismatches vs Stage B) — cross-checked cut-by-cut against the MAT/C++
source on 2026-06-12 (docs/open_questions.md "Reco-selection cross-check"):
CCInclusiveCuts.h + CVUniverse.h/MuonFunctions.h, instantiated with
minZ=5980, maxZ=8422, apothem=850, angle=20 deg, NoDeadtime(1) in
runEventLoop.cpp:360-368.

Selection (applied to data and MC reco trees identically):

  1. ZRange        5980 <= vtx[2] <= 8422 mm        (inclusive both ends)
  2. Apothem       (vtx[0], vtx[1]) inside the 850 mm flat-top hexagon (<)
  3. MaxMuonAngle  theta3d(muon_thetaX, muon_thetaY) < 20 deg (strict)
  4. HasMINOSMatch isMinosMatchTrack == 1
  5. NoDeadtime    phys_n_dead_discr_pair_upstream_prim_track_proj <= 1
  6. IsNeutrino    MasterAnaDev_minos_trk_qp < 0    (mu-)

Ordering note (load-bearing in the scalar chain, automatic here): the qp
charge sign is only physical for MINOS-matched tracks. In this vectorized
form the final selection is the AND of all masks, so cut 4 shields cut 6 by
construction; per-cut survivor counts (``cutflow``) apply the masks in the
certified order so they remain comparable to the scalar chain's cutflow.

All functions take NumPy arrays (uproot ``library="np"``) and return boolean
masks; no I/O. ``RECO_SELECTION_BRANCHES`` comes from config/branches.json
(role ``reco_selection``) — code cannot drift from the catalog.

Golden anchors (MC run 110040 / data run 10066, streamed gates in
tests/test_cuts.py): 844 selected data; 43643 selected MC, splitting into
43539 signal + 104 background (47 NC + 57 Wrong_Sign) with xsec.signal.
"""
import numpy as np

from . import catalog
from .constants import (DEAD_MAX_DISCR_PAIRS, MAX_MU_THETA_RAD, Z_MAX_MM,
                        Z_MIN_MM)
from .kinematics import in_hex_apothem, theta3d

RECO_SELECTION_BRANCHES = catalog.branches_for_role("reco_selection")

CUT_LABELS = (
    "ZRange",
    "Apothem",
    "MaxMuonAngle",
    "HasMINOSMatch",
    "NoDeadtime",
    "IsNeutrino",
)


def cut_z_range(vtx):
    """Tracker fiducial z, inclusive: 5980 <= z <= 8422 mm."""
    z = np.asarray(vtx, dtype=np.float64)[:, 2]
    return (Z_MIN_MM <= z) & (z <= Z_MAX_MM)


def cut_apothem(vtx):
    """Hexagonal transverse fiducial, strict <, apothem 850 mm."""
    v = np.asarray(vtx, dtype=np.float64)
    return in_hex_apothem(v[:, 0], v[:, 1])


def cut_max_muon_angle(muon_thetaX, muon_thetaY):
    """3D muon angle strictly below 20 degrees."""
    return theta3d(muon_thetaX, muon_thetaY) < MAX_MU_THETA_RAD


def cut_has_minos_match(isMinosMatchTrack):
    """Muon track matched to MINOS."""
    return np.asarray(isMinosMatchTrack) == 1


def cut_no_deadtime(n_dead_discr_pairs):
    """At most DEAD_MAX_DISCR_PAIRS dead discriminator pairs upstream."""
    return np.asarray(n_dead_discr_pairs) <= DEAD_MAX_DISCR_PAIRS


def cut_is_neutrino(minos_trk_qp):
    """Negative MINOS curvature (mu-). Only physical for matched tracks —
    always combine with cut_has_minos_match (reco_selection does)."""
    return np.asarray(minos_trk_qp) < 0


def _masks(arrays):
    """The six masks, keyed by CUT_LABELS, in the certified order."""
    return {
        "ZRange": cut_z_range(arrays["vtx"]),
        "Apothem": cut_apothem(arrays["vtx"]),
        "MaxMuonAngle": cut_max_muon_angle(arrays["muon_thetaX"],
                                           arrays["muon_thetaY"]),
        "HasMINOSMatch": cut_has_minos_match(arrays["isMinosMatchTrack"]),
        "NoDeadtime": cut_no_deadtime(
            arrays["phys_n_dead_discr_pair_upstream_prim_track_proj"]),
        "IsNeutrino": cut_is_neutrino(arrays["MasterAnaDev_minos_trk_qp"]),
    }


def reco_selection(arrays):
    """Final selection mask: AND of all six cuts.

    ``arrays`` is a mapping branch-name -> NumPy array containing (at least)
    RECO_SELECTION_BRANCHES, e.g. ``tree.arrays(RECO_SELECTION_BRANCHES,
    library="np")``.
    """
    masks = _masks(arrays)
    out = masks[CUT_LABELS[0]].copy()
    for label in CUT_LABELS[1:]:
        out &= masks[label]
    return out


def cutflow(arrays):
    """Survivor count after each cut, applied in the certified order.

    Returns a list of (label, n_surviving) pairs; the last count equals
    ``reco_selection(arrays).sum()``. Comparable to the scalar chain's
    cutflow because the order matches CUT_LABELS.
    """
    masks = _masks(arrays)
    flow, running = [], None
    for label in CUT_LABELS:
        running = masks[label] if running is None else (running & masks[label])
        flow.append((label, int(running.sum())))
    return flow
