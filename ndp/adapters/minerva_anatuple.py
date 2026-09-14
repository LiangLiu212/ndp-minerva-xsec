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

`build_cache` writes the per-file .npz tables every other stage reads (`ndp data cache`):
    truth_<tag>.npz             TruthTable of the Truth tree (MC only)
    reco_<tag>.npz              reco table: selection + RECO_CACHE_COLUMNS (data and MC)
    reco_<tag>_truthcols.npz    TruthTable of the reco rows' truth branches (MC only)
with <tag> = data10066 / mc110040 derived from the file name.
"""
from __future__ import annotations

import json
import math
import re
import time
from pathlib import Path

import numpy as np

from ..events import TruthTable, Normalization

MEV = 1e-3
MEV2 = 1e-6

#: bump when the cached reco columns change; readers refuse older caches.
RECO_CACHE_VERSION = 3

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
#: extra reco branches cached for user-defined observables: (branch, cache column, scale to GeV)
RECO_EXTRA_BRANCHES = (
    ("MasterAnaDev_minos_trk_p", "reco_minos_p", MEV),
    ("MasterAnaDev_recoil_E", "reco_recoil_E", MEV),                       # CC-inclusive recoil (== hadron_recoil_CCInc)
    ("MasterAnaDev_recoil_passivecorrected", "reco_recoil_E_passive", MEV),  # == hadron_recoil_default
    ("blob_recoil_E", "reco_recoil_E_calo", MEV),                          # == recoil_E_nopolyline
    ("recoil_E_polylinecorrected", "reco_recoil_E_polyline", MEV),
    ("MasterAnaDev_E", "reco_anatool_E_nu", MEV),
    ("MasterAnaDev_W", "reco_anatool_W", MEV),
    ("MasterAnaDev_x", "reco_anatool_x", 1.0),
    ("MasterAnaDev_y", "reco_anatool_y", 1.0),
    ("MasterAnaDev_visible_E", "reco_visible_E", MEV),
    ("n_prongs", "reco_n_prongs", 1.0),
    ("MasterAnaDev_hadron_number", "reco_n_hadron_tracks", 1.0),
    # --- cache v3: leading proton candidate, Michel electrons, isolated blobs (1mu1p selections) ---
    # MasterAnaDev_proton_* are event-level scalars of the tool's primary proton candidate (the
    # highest-momentum one: MasterAnaDev_sec_protons_P_fromdEdx never exceeds it on the open-data MC).
    # Momentum components are detector-frame (like MasterAnaDev_leptonE); MasterAnaDev_proton_theta is
    # beam-frame (== the rotated vector's angle, like muon_thetaX/Y). Sentinels (-9999 / -1) -> NaN.
    ("MasterAnaDev_proton_P_fromdEdx", "reco_proton_p", MEV),
    ("MasterAnaDev_proton_Px_fromdEdx", "reco_proton_px", MEV),
    ("MasterAnaDev_proton_Py_fromdEdx", "reco_proton_py", MEV),
    ("MasterAnaDev_proton_Pz_fromdEdx", "reco_proton_pz", MEV),
    ("MasterAnaDev_proton_E_fromdEdx", "reco_proton_E", MEV),
    ("MasterAnaDev_proton_T_fromdEdx", "reco_proton_T", MEV),
    ("MasterAnaDev_proton_theta", "reco_proton_theta_beam", 1.0),
    ("MasterAnaDev_proton_score1", "reco_proton_score1", 1.0),
    ("MasterAnaDev_proton_score2", "reco_proton_score2", 1.0),
    ("MasterAnaDev_pion_score1", "reco_pion_score1", 1.0),
    ("MasterAnaDev_proton_startPointZ", "reco_proton_start_z", 1.0),
    ("MasterAnaDev_proton_endPointZ", "reco_proton_end_z", 1.0),
    ("MasterAnaDev_proton_patternRec", "reco_proton_pattern_rec", 1.0),
    ("MasterAnaDev_sec_protons_P_fromdEdx_sz", "reco_n_sec_protons", 1.0),
    ("improved_nmichel", "reco_n_michel", 1.0),
    ("n_nonvtx_iso_blobs", "reco_n_iso_blobs", 1.0),
    ("n_nonvtx_iso_blobs_all", "reco_n_iso_blobs_all", 1.0),
    ("nonvtx_iso_blobs_energy", "reco_iso_blobs_E", MEV),
    ("proton_prong_PDG", "reco_proton_truth_pdg", 1.0),            # MC: truth PDG of the candidate; data: -1
)
#: columns whose tuple sentinel (-9999, or -1 for the end point) means "no candidate" -> NaN in the cache
RECO_SENTINEL_NAN = ("reco_proton_p", "reco_proton_px", "reco_proton_py", "reco_proton_pz", "reco_proton_E", "reco_proton_T",
                     "reco_proton_theta_beam", "reco_proton_score1", "reco_proton_score2", "reco_pion_score1",
                     "reco_proton_start_z", "reco_proton_end_z")
#: fixed-size array branches cached as one element: (branch, index, cache column)
RECO_ARRAY_ELEMENTS = (("MasterAnaDev_hadron_isExiting", 0, "reco_proton_exiting"),)
RECO_CACHE_COLUMNS = ("passed", "failing_cut", "reco_minos_matched", "reco_p", "reco_theta", "reco_pT", "reco_pz", "reco_E_mu",
                      "reco_thetaX", "reco_thetaY", "reco_minos_qp", "reco_vtx_x", "reco_vtx_y", "reco_vtx_z",
                      "reco_mu_px", "reco_mu_py", "reco_mu_pz", "reco_n_dead_discr",
                      *[c for _, c, _ in RECO_EXTRA_BRANCHES], *[c for _, _, c in RECO_ARRAY_ELEMENTS])
TRUTH_SCALARS = ("mc_incoming", "mc_current", "mc_intType", "mc_targetZ", "mc_targetA",
                 "mc_incomingE", "mc_Q2", "mc_w", "mc_primaryLepton")
TRUTH_VECTORS = ("mc_primFSLepton", "mc_vtx")
TRUTH_FS = ("mc_nFSPart", "mc_FSPartPDG", "mc_FSPartE", "mc_FSPartPx", "mc_FSPartPy", "mc_FSPartPz")

#: MINERvA `mc_intType` -> NDP interaction code. MINERvA re-labels GENIE's enum:
#: 1=QE, 2=RES, 3=DIS, 4=COH, 8=MEC (2p2h); 5-7 are rare electron-scattering / IMD types.
MINERVA_INT_TYPE = {1: 1, 2: 2, 3: 3, 4: 4, 8: 5}

_TAG = re.compile(r"_(data|mc)_AnaTuple_run0*(\d+)")


def cache_tag(filename: str | Path) -> str:
    """MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root -> 'mc110040' (data -> 'data10066')."""
    m = _TAG.search(Path(filename).name)
    if not m:
        raise ValueError(f"cannot derive a cache tag from {filename!r}")
    return f"{m.group(1)}{m.group(2)}"


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


# --- CCQE-like 1mu1p selection (arXiv:2503.15047 Sec. "Analysis and results"), on cached columns ---
CCQELIKE_1MU1P_LABELS = ("ZRange", "Apothem", "HasMINOSMatch", "NoDeadtime", "IsNeutrino", "MuonWindow",
                         "HasProtonCandidate", "ProtonContained", "ProtonScore", "ProtonWindow", "NoMichel", "IsoBlobs")


def ccqelike_1mu1p_cutflow(r: dict, params: dict) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised 1mu1p chain on a v3 reco cache. Returns (passed, failing_cut_index; -1 = passed).

    params (all from the channel manifest, none defaulted here):
        muon:   {theta_max_deg, p_min_gev, p_max_gev}          reco muon window (beam-frame angle)
        proton: {theta_max_deg, p_min_gev, p_max_gev}          reco proton window (beam-frame angle, dE/dx momentum)
        proton_score1_min, require_contained (bool), n_michel_max, n_iso_blobs_max
    A NaN kinematic (no candidate) fails the cut it is tested in, as in the tool.
    """
    mu, pr = params["muon"], params["proton"]
    vx, vy, vz = r["reco_vtx_x"], r["reco_vtx_y"], r["reco_vtx_z"]
    th_mu, p_mu = np.asarray(r["reco_theta"], float), np.asarray(r["reco_p"], float)
    p_p, th_p = np.asarray(r["reco_proton_p"], float), np.asarray(r["reco_proton_theta_beam"], float)
    score = np.asarray(r["reco_proton_score1"], float)
    with np.errstate(invalid="ignore"):
        fail = [
            ~((Z_MIN <= vz) & (vz <= Z_MAX)),
            ~((np.abs(vx) < APOTHEM) & (np.abs(vy) < APOTHEM_SLOPE * np.abs(vx) + APOTHEM_INTERCEPT)),
            ~np.asarray(r["reco_minos_matched"], bool),
            np.asarray(r["reco_n_dead_discr"], float) > DEAD_MAX,
            ~(np.asarray(r["reco_minos_qp"], float) < 0),
            ~((th_mu < np.deg2rad(mu["theta_max_deg"])) & (p_mu > mu["p_min_gev"]) & (p_mu < mu["p_max_gev"])),
            ~(p_p > 0),
            (np.asarray(r["reco_proton_exiting"], float) == 1) if params.get("require_contained", True) else np.zeros(len(vz), bool),
            ~(score > params["proton_score1_min"]),
            ~((th_p < np.deg2rad(pr["theta_max_deg"])) & (p_p > pr["p_min_gev"]) & (p_p < pr["p_max_gev"])),
            ~(np.asarray(r["reco_n_michel"], float) <= params["n_michel_max"]),
            ~(np.asarray(r["reco_n_iso_blobs"], float) <= params["n_iso_blobs_max"]),
        ]
    failing = np.full(len(vz), -1, dtype=np.int64)
    alive = np.ones(len(vz), dtype=bool)
    for i, f in enumerate(fail):
        killed = alive & f
        failing[killed] = i
        alive &= ~f
    return alive, failing


def reco_muon_kinematics(a: dict) -> dict:
    """Reco muon (p, theta, pT, pz, E) in GeV from MasterAnaDev_leptonE + projected angles."""
    le = a["MasterAnaDev_leptonE"]
    p = np.sqrt(le[:, 0] ** 2 + le[:, 1] ** 2 + le[:, 2] ** 2) * MEV
    th = theta3d(a["muon_thetaX"], a["muon_thetaY"])
    return {"p": p, "theta": th, "pT": p * np.sin(th), "pz": p * np.cos(th), "E": le[:, 3] * MEV}


def reco_columns(a: dict) -> dict:
    """The full cached reco table (RECO_CACHE_COLUMNS) from the raw branch arrays."""
    passed, failing = cc_inclusive_cutflow(a)
    k = reco_muon_kinematics(a)
    vtx = a["vtx"]
    out = {"passed": passed, "failing_cut": failing, "reco_minos_matched": np.asarray(a["isMinosMatchTrack"]) == 1,
           "reco_p": k["p"], "reco_theta": k["theta"], "reco_pT": k["pT"],
           "reco_pz": k["pz"], "reco_E_mu": k["E"], "reco_thetaX": np.asarray(a["muon_thetaX"], float),
           "reco_thetaY": np.asarray(a["muon_thetaY"], float), "reco_minos_qp": np.asarray(a["MasterAnaDev_minos_trk_qp"], float),
           "reco_vtx_x": vtx[:, 0].astype(float), "reco_vtx_y": vtx[:, 1].astype(float), "reco_vtx_z": vtx[:, 2].astype(float)}
    le = a["MasterAnaDev_leptonE"]
    out["reco_mu_px"], out["reco_mu_py"], out["reco_mu_pz"] = le[:, 0] * MEV, le[:, 1] * MEV, le[:, 2] * MEV   # detector frame
    out["reco_n_dead_discr"] = np.asarray(a["phys_n_dead_discr_pair_upstream_prim_track_proj"], float)
    for branch, col, scale in RECO_EXTRA_BRANCHES:
        if branch in a:
            v = np.asarray(a[branch], float)
            if col in RECO_SENTINEL_NAN:
                v = np.where(v <= -999.0, np.nan, v)
                if col == "reco_proton_end_z":
                    v = np.where(v == -1.0, np.nan, v)
            out[col] = v * scale
    for branch, idx, col in RECO_ARRAY_ELEMENTS:
        if branch in a:
            arr = np.asarray(a[branch])
            out[col] = arr[:, idx].astype(float) if arr.ndim == 2 else np.full(len(vtx), np.nan)
    return out


def read_reco(path: str | Path, is_mc: bool, entry_stop: int | None = None, extra: bool = True, with_fs: bool = True) -> dict:
    """Read the `MasterAnaDev` tree -> dict of arrays with selection + kinematics.

    Keys: `passed`, `failing_cut`, `reco` (dict p/theta/pT/pz/E), `columns` (the cache
    table, see RECO_CACHE_COLUMNS), and for MC `truth` (a TruthTable built from the reco
    rows' truth branches, in GeV).
    """
    tree = _uproot().open(path)["MasterAnaDev"]
    keys = set(tree.keys())
    branches = list(RECO_BRANCHES)
    if extra:
        branches += [b for b, _, _ in RECO_EXTRA_BRANCHES if b in keys]
        branches += [b for b, _, _ in RECO_ARRAY_ELEMENTS if b in keys]
    if is_mc:
        branches += list(TRUTH_SCALARS) + list(TRUTH_VECTORS)
    a = tree.arrays(branches, library="np", entry_stop=entry_stop)
    cols = reco_columns(a)
    out = {"n_entries": len(cols["passed"]), "passed": cols["passed"], "failing_cut": cols["failing_cut"],
           "reco": reco_muon_kinematics(a), "columns": cols, "path": str(path),
           "missing_extra_branches": [b for b, _, _ in RECO_EXTRA_BRANCHES if b not in keys]}
    if is_mc:
        fs = _read_fs(tree, entry_stop) if with_fs else None
        out["truth"] = _truth_table_from_arrays(a, source=f"{path}:MasterAnaDev", fs=fs)
    return out


def load_reco_cache(path: str | Path, require_version: int | None = RECO_CACHE_VERSION) -> dict:
    """reco_<tag>.npz -> dict of arrays (+ '__meta__' dict). Older caches raise with the fix."""
    z = np.load(path, allow_pickle=False)
    r = {k: z[k] for k in z.files if k != "__meta__"}
    meta = json.loads(str(z["__meta__"])) if "__meta__" in z.files else {}
    v = int(meta.get("cache_version", 1))
    if require_version is not None and v < require_version:
        raise ValueError(f"{path} is reco-cache version {v} (< {require_version}): rebuild it with `python -m ndp data cache`")
    r["__meta__"] = meta
    return r


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


def _read_fs(tree, entry_stop: int | None = None) -> dict:
    """Jagged mc_FSPart* branches of a tree -> CSR final-state columns in GeV."""
    import awkward as ak
    parts = tree.arrays(["mc_FSPartPDG", "mc_FSPartE", "mc_FSPartPx", "mc_FSPartPy", "mc_FSPartPz"],
                        library="ak", entry_stop=entry_stop)
    counts = ak.num(parts["mc_FSPartPDG"], axis=1).to_numpy()
    return {
        "fs_offsets": np.concatenate([[0], np.cumsum(counts)]).astype(np.int64),
        "fs_pdg": ak.flatten(parts["mc_FSPartPDG"]).to_numpy(),
        "fs_E": ak.flatten(parts["mc_FSPartE"]).to_numpy() * MEV,
        "fs_px": ak.flatten(parts["mc_FSPartPx"]).to_numpy() * MEV,
        "fs_py": ak.flatten(parts["mc_FSPartPy"]).to_numpy() * MEV,
        "fs_pz": ak.flatten(parts["mc_FSPartPz"]).to_numpy() * MEV,
    }


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
    fs = _read_fs(tree, entry_stop) if with_fs else None
    pot = read_pot(path)["pot_used"]
    t = _truth_table_from_arrays(a, source=f"{path}:Truth", fs=fs, pot=pot)
    t.meta["n_truth_entries_in_file"] = int(tree.num_entries)
    t.meta["entry_stop"] = entry_stop
    return t


# --------------------------------------------------------------------------------------
# Cache builder (`ndp data cache`)
# --------------------------------------------------------------------------------------
def build_cache(path: str | Path, cache_dir: str | Path, is_mc: bool, *, truth: bool = True, reco: bool = True,
                entry_stop: int | None = None, log=print) -> dict:
    """Write the .npz tables for one AnaTuple. Returns what was written with timings."""
    from ..io import cheap_fingerprint
    path, cache_dir = Path(path), Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    tag = cache_tag(path)
    out: dict = {"tag": tag, "path": str(path), "written": []}
    if is_mc and truth:
        t0 = time.time()
        tt = read_truth(path, entry_stop=entry_stop)
        p = cache_dir / f"truth_{tag}.npz"
        tt.save(p); out["written"].append(str(p)); out["truth_s"] = round(time.time() - t0, 1); out["n_truth"] = tt.n
        log(f"truth {tag}: {tt.n} entries in {out['truth_s']} s -> {p}")
    if reco:
        t0 = time.time()
        r = read_reco(path, is_mc=is_mc, entry_stop=entry_stop)
        cols = r["columns"]
        meta = {"cache_version": RECO_CACHE_VERSION, "source": str(path), "tree": "MasterAnaDev", "is_mc": is_mc,
                "units": "GeV, GeV^2, mm", "selection": "minerva_cc_inclusive_v1 (vectorised cc_inclusive_cutflow)",
                "columns": list(cols), "branches": {c: b for b, c, _ in RECO_EXTRA_BRANCHES if c in cols},
                "missing_extra_branches": r["missing_extra_branches"], "entry_stop": entry_stop,
                "n_entries": int(r["n_entries"]), "n_passed": int(cols["passed"].sum()),
                "file_fingerprint": cheap_fingerprint(path), "built": time.strftime("%Y-%m-%dT%H:%M:%S")}
        p = cache_dir / f"reco_{tag}.npz"
        np.savez_compressed(p, __meta__=json.dumps(meta), **cols)
        out["written"].append(str(p))
        if is_mc:
            p2 = cache_dir / f"reco_{tag}_truthcols.npz"
            r["truth"].save(p2); out["written"].append(str(p2))
        out["reco_s"] = round(time.time() - t0, 1); out["n_reco"] = int(r["n_entries"]); out["n_passed"] = int(cols["passed"].sum())
        log(f"reco {tag}: {r['n_entries']} entries, {out['n_passed']} selected, in {out['reco_s']} s -> {p}")
    return out


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
