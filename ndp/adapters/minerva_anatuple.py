"""MINERvA MasterAnaDev AnaTuple adapter (uproot, vectorised).

Two trees matter:
  * `Truth`        — every generated interaction in the simulated volume (the efficiency
                     denominator and the reference generator sample; GENIE truth only).
  * `MasterAnaDev` — one row per reconstructed candidate, with reco branches (data + MC)
                     and, for MC, the matching truth branches.

Units in the tuple are MeV / MeV^2 / mm; this adapter converts to GeV / GeV^2 / mm.

The reconstruction-level selection here is a vectorised transcription of
`tools/cc_inclusive_selector.py` from the MINERvA exploration repo, cut for cut and
NaN behaviour included; `parity_vs_tool` certifies the two agree row by row.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from ..events import TruthTable, Normalization

MEV = 1e-3
MEV2 = 1e-6

# --- selection constants (MINERvA-101 tutorial; mirror of tools/cc_inclusive_selector.py) ----
Z_MIN, Z_MAX = 5980.0, 8422.0
APOTHEM = 850.0
APOTHEM_SLOPE = -1.0 / math.sqrt(3.0)
APOTHEM_INTERCEPT = 2.0 * APOTHEM / math.sqrt(3.0)
MAX_MU_THETA_RAD = 20.0 * math.pi / 180.0
DEAD_MAX = 1
CUT_LABELS = ("ZRange", "Apothem", "MaxMuonAngle", "HasMINOSMatch", "NoDeadtime", "IsNeutrino")

RECO_BRANCHES = ("vtx", "muon_thetaX", "muon_thetaY", "isMinosMatchTrack",
                 "phys_n_dead_discr_pair_upstream_prim_track_proj", "MasterAnaDev_minos_trk_qp",
                 "MasterAnaDev_leptonE")
TRUTH_SCALARS = ("mc_incoming", "mc_current", "mc_intType", "mc_targetZ", "mc_targetA",
                 "mc_incomingE", "mc_Q2", "mc_w", "mc_primaryLepton")
TRUTH_VECTORS = ("mc_primFSLepton", "mc_vtx")
TRUTH_FS = ("mc_nFSPart", "mc_FSPartPDG", "mc_FSPartE", "mc_FSPartPx", "mc_FSPartPy", "mc_FSPartPz")

#: MINERvA `mc_intType` -> NDP interaction code. MINERvA re-labels GENIE's enum:
#: 1=QE, 2=RES, 3=DIS, 4=COH, 8=MEC (2p2h); 5-7 are rare electron-scattering / IMD types.
MINERVA_INT_TYPE = {1: 1, 2: 2, 3: 3, 4: 4, 8: 5}


def _uproot():
    import uproot  # local import so the rest of the package imports without uproot
    return uproot


def theta3d(theta_x: np.ndarray, theta_y: np.ndarray) -> np.ndarray:
    tx, ty = np.tan(theta_x), np.tan(theta_y)
    return np.arccos(1.0 / np.sqrt(1.0 + tx * tx + ty * ty))


def read_pot(path: str | Path) -> dict:
    t = _uproot().open(path)["Meta"]
    a = t.arrays(["POT_Used", "POT_Total"], library="np")
    return {"pot_used": float(np.sum(a["POT_Used"])), "pot_total": float(np.sum(a["POT_Total"])),
            "n_meta_entries": int(len(a["POT_Used"]))}


# --------------------------------------------------------------------------------------
# Reconstruction level
# --------------------------------------------------------------------------------------
def cc_inclusive_cutflow(a: dict) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised 6-cut chain. Returns (passed[bool], failing_cut_index[int], -1 = passed).

    Each `fail_*` expression is the literal negation-free transcription of the tool's
    `if not (...)` / `if (...)` test so NaN semantics match (a NaN muon angle passes the
    angle cut in both implementations because `nan >= x` is False).
    """
    vtx = a["vtx"]
    vx, vy, vz = vtx[:, 0], vtx[:, 1], vtx[:, 2]
    fail = [
        ~((Z_MIN <= vz) & (vz <= Z_MAX)),
        ~((np.abs(vx) < APOTHEM) & (np.abs(vy) < APOTHEM_SLOPE * np.abs(vx) + APOTHEM_INTERCEPT)),
        theta3d(a["muon_thetaX"], a["muon_thetaY"]) >= MAX_MU_THETA_RAD,
        a["isMinosMatchTrack"] != 1,
        a["phys_n_dead_discr_pair_upstream_prim_track_proj"] > DEAD_MAX,
        ~(a["MasterAnaDev_minos_trk_qp"] < 0),
    ]
    failing = np.full(len(vz), -1, dtype=np.int64)
    alive = np.ones(len(vz), dtype=bool)
    for i, f in enumerate(fail):
        killed = alive & f
        failing[killed] = i
        alive &= ~f
    return alive, failing


def reco_muon_kinematics(a: dict) -> dict:
    """Reco muon (p, theta, pT, pz) in GeV from MasterAnaDev_leptonE + projected angles."""
    le = a["MasterAnaDev_leptonE"]
    p = np.sqrt(le[:, 0] ** 2 + le[:, 1] ** 2 + le[:, 2] ** 2) * MEV
    th = theta3d(a["muon_thetaX"], a["muon_thetaY"])
    return {"p": p, "theta": th, "pT": p * np.sin(th), "pz": p * np.cos(th)}


def read_reco(path: str | Path, is_mc: bool, entry_stop: int | None = None) -> dict:
    """Read the `MasterAnaDev` tree -> dict of arrays with selection + kinematics.

    Keys: `passed`, `failing_cut`, `reco` (dict p/theta/pT/pz), and for MC `truth`
    (a TruthTable built from the reco rows' truth branches, in GeV).
    """
    tree = _uproot().open(path)["MasterAnaDev"]
    branches = list(RECO_BRANCHES)
    if is_mc:
        branches += list(TRUTH_SCALARS) + list(TRUTH_VECTORS)
    a = tree.arrays(branches, library="np", entry_stop=entry_stop)
    passed, failing = cc_inclusive_cutflow(a)
    out = {"n_entries": len(passed), "passed": passed, "failing_cut": failing,
           "reco": reco_muon_kinematics(a), "path": str(path)}
    if is_mc:
        out["truth"] = _truth_table_from_arrays(a, source=f"{path}:MasterAnaDev")
    return out


# --------------------------------------------------------------------------------------
# Truth level
# --------------------------------------------------------------------------------------
def _truth_table_from_arrays(a: dict, source: str, fs: dict | None = None,
                             pot: float | None = None) -> TruthTable:
    fl = a["mc_primFSLepton"]
    v = a["mc_vtx"]
    it = a["mc_intType"]
    int_type = np.zeros(len(it), dtype=np.int64)
    for k, code in MINERVA_INT_TYPE.items():
        int_type[it == k] = code
    cols = {
        "nu_pdg": a["mc_incoming"], "E_nu": a["mc_incomingE"] * MEV,
        "lep_pdg": a.get("mc_primaryLepton", np.where(a["mc_current"] == 1, np.sign(a["mc_incoming"]) * 13, a["mc_incoming"])),
        "lep_px": fl[:, 0] * MEV, "lep_py": fl[:, 1] * MEV, "lep_pz": fl[:, 2] * MEV, "lep_E": fl[:, 3] * MEV,
        "current": a["mc_current"], "int_type": int_type,
        "target_Z": a["mc_targetZ"], "target_A": a["mc_targetA"],
        "Q2": a["mc_Q2"] * MEV2, "W": a["mc_w"] * MEV,
        "weight": np.ones(len(it)),
        "vtx_x": v[:, 0], "vtx_y": v[:, 1], "vtx_z": v[:, 2],
        "generator_int_type": it,
    }
    if fs is not None:
        cols.update(fs)
    meta = {"source": source, "generator": "GENIE 2.12.6 (MINERvA Open Data StandardMC; unconfirmed tag)",
            "units": "GeV, GeV^2, mm", "frame": "detector", "has_geometry": True,
            "norm": Normalization(kind="pot", pot=pot, notes="MINERvA MC exposure from the Meta tree").to_dict()
            if pot else Normalization(kind="shape").to_dict()}
    return TruthTable(cols, meta)


def read_truth(path: str | Path, entry_stop: int | None = None, with_fs: bool = True,
               step_size: str = "200 MB") -> TruthTable:
    """Read the `Truth` tree (all generated events) into a TruthTable, GeV units.

    Final-state particles are read jagged and stored CSR-style. ~90 s for the full
    544k-entry me1A MC file on the EAF filesystem; callers should cache the .npz.
    """
    up = _uproot()
    f = up.open(path)
    tree = f["Truth"]
    branches = list(TRUTH_SCALARS) + list(TRUTH_VECTORS)
    a = tree.arrays(branches, library="np", entry_stop=entry_stop)
    fs = None
    if with_fs:
        import awkward as ak
        parts = tree.arrays(["mc_FSPartPDG", "mc_FSPartE", "mc_FSPartPx", "mc_FSPartPy", "mc_FSPartPz"],
                            library="ak", entry_stop=entry_stop)
        counts = ak.num(parts["mc_FSPartPDG"], axis=1).to_numpy()
        fs = {
            "fs_offsets": np.concatenate([[0], np.cumsum(counts)]).astype(np.int64),
            "fs_pdg": ak.flatten(parts["mc_FSPartPDG"]).to_numpy(),
            "fs_E": ak.flatten(parts["mc_FSPartE"]).to_numpy() * MEV,
            "fs_px": ak.flatten(parts["mc_FSPartPx"]).to_numpy() * MEV,
            "fs_py": ak.flatten(parts["mc_FSPartPy"]).to_numpy() * MEV,
            "fs_pz": ak.flatten(parts["mc_FSPartPz"]).to_numpy() * MEV,
        }
    pot = read_pot(path)["pot_used"]
    t = _truth_table_from_arrays(a, source=f"{path}:Truth", fs=fs, pot=pot)
    t.meta["n_truth_entries_in_file"] = int(tree.num_entries)
    t.meta["entry_stop"] = entry_stop
    return t


# --------------------------------------------------------------------------------------
# Certification helper
# --------------------------------------------------------------------------------------
def parity_vs_tool(path: str | Path, minerva_repo: str | Path, n: int = 20000) -> dict:
    """Compare this module's vectorised selection with the MINERvA repo's per-entry tool.

    Builds one lightweight row object per entry (the tool reads attributes named after
    branches) and checks pass/fail AND the failing-cut label agree for every row.
    """
    import sys
    sys.path.insert(0, str(minerva_repo))
    from tools.cc_inclusive_selector import passes_cc_inclusive, CUT_LABELS as TOOL_LABELS  # type: ignore

    tree = _uproot().open(path)["MasterAnaDev"]
    a = tree.arrays(list(RECO_BRANCHES[:-1]), library="np", entry_stop=n)
    passed, failing = cc_inclusive_cutflow(a)

    class Row:  # the tool only needs attribute access
        __slots__ = ("vtx", "muon_thetaX", "muon_thetaY", "isMinosMatchTrack",
                     "phys_n_dead_discr_pair_upstream_prim_track_proj", "MasterAnaDev_minos_trk_qp")

    mism = []
    for i in range(len(passed)):
        r = Row()
        r.vtx = a["vtx"][i]
        r.muon_thetaX = float(a["muon_thetaX"][i]); r.muon_thetaY = float(a["muon_thetaY"][i])
        r.isMinosMatchTrack = int(a["isMinosMatchTrack"][i])
        r.phys_n_dead_discr_pair_upstream_prim_track_proj = int(a["phys_n_dead_discr_pair_upstream_prim_track_proj"][i])
        r.MasterAnaDev_minos_trk_qp = float(a["MasterAnaDev_minos_trk_qp"][i])
        ok, cut = passes_cc_inclusive(r)
        mine = None if failing[i] < 0 else CUT_LABELS[failing[i]]
        if ok != bool(passed[i]) or cut != mine:
            mism.append((i, ok, cut, bool(passed[i]), mine))
    return {"n_compared": int(len(passed)), "n_mismatch": len(mism), "first_mismatches": mism[:5],
            "tool_labels": list(TOOL_LABELS), "n_passed": int(passed.sum())}
