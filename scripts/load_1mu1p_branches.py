#!/usr/bin/env python3
"""Load the AnaTuple branches used by the 1mu + leading proton analysis with uproot.

Standalone (numpy + awkward + uproot only). The branch lists mirror
`ndp/adapters/minerva_anatuple.py` (RECO_BRANCHES, RECO_EXTRA_BRANCHES, RECO_ARRAY_ELEMENTS, TRUTH_*).
Values are returned in the tuple's own units (MeV, mm, radians); nothing is rescaled or NaN-ified.

Usage (repo pixi env: uproot + XRootD):
    pixi run python scripts/load_1mu1p_branches.py <file.root | root://...> [--entry-stop N] [--out file.npz]

From Python:
    from load_1mu1p_branches import load_anatuple
    ev = load_anatuple("MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root", entry_stop=10000)
    ev["reco"]["MasterAnaDev_proton_P_fromdEdx"]      # numpy array, one entry per reco candidate
    ev["truth"]["mc_FSPartPDG"]                        # awkward jagged array, one entry per true event (MC)
    ev["pot"]["pot_used"]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# --- MasterAnaDev tree: reco branches (data and MC) ---------------------------------------------
#: used directly by the selection cuts
RECO_SELECTION = (
    "vtx",                                                  # [x, y, z] mm; fiducial z range + hexagon apothem
    "muon_thetaX", "muon_thetaY",                           # beam-frame muon angles (rad)
    "isMinosMatchTrack",
    "phys_n_dead_discr_pair_upstream_prim_track_proj",      # dead time
    "MasterAnaDev_minos_trk_qp",                            # q/p in MINOS: < 0 selects mu-
    "MasterAnaDev_leptonE",                                 # [px, py, pz, E] MeV, detector frame
    "MasterAnaDev_proton_P_fromdEdx",                       # MeV/c, leading proton candidate
    "MasterAnaDev_proton_theta",                            # beam-frame angle (rad)
    "MasterAnaDev_proton_score1",
    "improved_nmichel",
    "n_nonvtx_iso_blobs",
)
#: cached for observables and diagnostics
RECO_EXTRA = (
    "MasterAnaDev_proton_Px_fromdEdx", "MasterAnaDev_proton_Py_fromdEdx", "MasterAnaDev_proton_Pz_fromdEdx",
    "MasterAnaDev_proton_E_fromdEdx", "MasterAnaDev_proton_T_fromdEdx",
    "MasterAnaDev_proton_score2", "MasterAnaDev_pion_score1",
    "MasterAnaDev_proton_startPointZ", "MasterAnaDev_proton_endPointZ", "MasterAnaDev_proton_patternRec",
    "MasterAnaDev_sec_protons_P_fromdEdx_sz",
    "proton_prong_PDG",                                     # MC: truth PDG of the candidate; data: -1
    "n_nonvtx_iso_blobs_all", "nonvtx_iso_blobs_energy", "n_prongs", "MasterAnaDev_hadron_number",
    "MasterAnaDev_minos_trk_p",
    "MasterAnaDev_recoil_E", "MasterAnaDev_recoil_passivecorrected", "blob_recoil_E", "recoil_E_polylinecorrected",
    "MasterAnaDev_visible_E",
    "MasterAnaDev_E", "MasterAnaDev_W", "MasterAnaDev_x", "MasterAnaDev_y",
)
#: fixed-size array branches of which one element is used: (branch, index)
RECO_ARRAY_ELEMENTS = (("MasterAnaDev_hadron_isExiting", 0),)

# --- truth branches (MC only; present in both the MasterAnaDev and the Truth tree) ---------------
TRUTH_SCALARS = ("mc_incoming", "mc_incomingE", "mc_current", "mc_intType", "mc_targetZ", "mc_targetA",
                 "mc_Q2", "mc_w", "mc_primaryLepton")
TRUTH_VECTORS = ("mc_primFSLepton", "mc_vtx")               # [px,py,pz,E] MeV and [x,y,z,t]
TRUTH_FS = ("mc_nFSPart", "mc_FSPartPDG", "mc_FSPartE", "mc_FSPartPx", "mc_FSPartPy", "mc_FSPartPz")

# --- Meta tree ---------------------------------------------------------------------------------
META_BRANCHES = ("POT_Used", "POT_Total", "Total_Reco_Entries", "Total_Truth_Entries")

_TAG = re.compile(r"_(data|mc)_AnaTuple_run0*(\d+)")


def is_mc_file(path: str) -> bool | None:
    """True/False from the AnaTuple file name, None if the name does not follow the open-data pattern."""
    m = _TAG.search(Path(str(path)).name)
    return None if m is None else m.group(1) == "mc"


def _present(tree, wanted) -> tuple[list[str], list[str]]:
    keys = set(tree.keys())
    have = [b for b in wanted if b in keys]
    return have, [b for b in wanted if b not in keys]


def read_meta(f) -> dict:
    """POT and entry totals summed over the Meta tree."""
    tree = f["Meta"]
    have, _ = _present(tree, META_BRANCHES)
    a = tree.arrays(have, library="np")
    out = {"pot_used": float(np.sum(a["POT_Used"])), "pot_total": float(np.sum(a["POT_Total"])),
           "n_meta_entries": int(len(a["POT_Used"]))}
    for k in ("Total_Reco_Entries", "Total_Truth_Entries"):
        if k in a:
            out[k.lower()] = int(np.sum(a[k]))
    return out


def read_reco(f, mc: bool, entry_stop: int | None = None) -> tuple[dict, list[str]]:
    """MasterAnaDev reco branches (+ the matching truth branches for MC) as numpy arrays."""
    tree = f["MasterAnaDev"]
    wanted = list(RECO_SELECTION) + list(RECO_EXTRA) + [b for b, _ in RECO_ARRAY_ELEMENTS]
    if mc:
        wanted += list(TRUTH_SCALARS) + list(TRUTH_VECTORS)
    have, missing = _present(tree, wanted)
    a = tree.arrays(have, library="np", entry_stop=entry_stop)
    for branch, idx in RECO_ARRAY_ELEMENTS:                 # keep the one element the analysis uses
        if branch in a:
            arr = a[branch]
            a[f"{branch}[{idx}]"] = arr[:, idx] if arr.ndim == 2 else np.array([row[idx] for row in arr])
    return a, missing


def read_truth(f, entry_stop: int | None = None) -> tuple[dict, list[str]]:
    """Truth tree (every generated event, MC only): scalars/vectors as numpy, mc_FSPart* as awkward jagged."""
    import awkward as ak
    tree = f["Truth"]
    have, missing = _present(tree, list(TRUTH_SCALARS) + list(TRUTH_VECTORS) + list(TRUTH_FS))
    flat = [b for b in have if b not in TRUTH_FS]
    jag = [b for b in have if b in TRUTH_FS]
    out = dict(tree.arrays(flat, library="np", entry_stop=entry_stop))
    if jag:
        parts = tree.arrays(jag, library="ak", entry_stop=entry_stop)
        for b in jag:
            out[b] = parts[b]                               # ak.Array; ak.flatten / ak.num for CSR form
        if "mc_FSPartPDG" in out:
            out["fs_counts"] = ak.num(out["mc_FSPartPDG"], axis=1).to_numpy()
    return out, missing


def load_anatuple(path: str, entry_stop: int | None = None, mc: bool | None = None, timeout: int = 300) -> dict:
    """Open a local path or root:// URL and load every branch of the 1mu1p analysis."""
    import uproot
    if mc is None:
        mc = is_mc_file(path)
    f = uproot.open(str(path), timeout=int(timeout))        # XRootD wants an integer timeout
    if mc is None:                                          # unknown name: look for truth in the tuple
        mc = "mc_incoming" in f["MasterAnaDev"].keys()
    ev = {"path": str(path), "is_mc": bool(mc), "pot": read_meta(f)}
    ev["reco"], ev["missing_reco"] = read_reco(f, mc, entry_stop)
    if mc and "Truth" in f:
        ev["truth"], ev["missing_truth"] = read_truth(f, entry_stop)
    else:
        ev["truth"], ev["missing_truth"] = {}, []
    return ev


def save_npz(ev: dict, out: str) -> None:
    """Flat npz: reco_<branch>, truth_<branch>; jagged FS branches become truth_<branch>_flat + truth_fs_counts."""
    import awkward as ak
    arrays = {}
    for k, v in ev["reco"].items():
        arrays[f"reco_{k}"] = v
    for k, v in ev["truth"].items():
        if isinstance(v, ak.Array):
            arrays[f"truth_{k}_flat"] = ak.flatten(v).to_numpy()
        else:
            arrays[f"truth_{k}"] = v
    arrays["pot_used"] = np.array(ev["pot"]["pot_used"])
    arrays["pot_total"] = np.array(ev["pot"]["pot_total"])
    np.savez(out, **arrays)


def _summary(ev: dict) -> str:
    lines = [f"{ev['path']}", f"  is_mc: {ev['is_mc']}", f"  POT: {ev['pot']}"]
    n = len(next(iter(ev["reco"].values()))) if ev["reco"] else 0
    lines.append(f"  MasterAnaDev: {n} rows, {len(ev['reco'])} arrays"
                 + (f", missing {ev['missing_reco']}" if ev["missing_reco"] else ""))
    for k in ("MasterAnaDev_proton_P_fromdEdx", "MasterAnaDev_leptonE", "vtx", "MasterAnaDev_hadron_isExiting[0]"):
        if k in ev["reco"]:
            a = ev["reco"][k]
            lines.append(f"    {k}: shape {getattr(a, 'shape', None)} dtype {a.dtype}")
    if ev["truth"]:
        nt = len(ev["truth"]["mc_incoming"]) if "mc_incoming" in ev["truth"] else "?"
        lines.append(f"  Truth: {nt} events, {len(ev['truth'])} arrays"
                     + (f", missing {ev['missing_truth']}" if ev["missing_truth"] else ""))
        if "fs_counts" in ev["truth"]:
            c = ev["truth"]["fs_counts"]
            lines.append(f"    mc_FSPart*: {int(c.sum())} particles, mean {c.mean():.2f} per event")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("path", help="AnaTuple file: local path or root://fndcadoor.fnal.gov:1095/... URL")
    ap.add_argument("--entry-stop", type=int, default=None, help="read only the first N entries of each tree")
    ap.add_argument("--mc", action="store_true", help="force MC handling (default: from the file name)")
    ap.add_argument("--data", action="store_true", help="force data handling")
    ap.add_argument("--out", default=None, help="save the arrays to this .npz")
    args = ap.parse_args(argv)
    mc = True if args.mc else (False if args.data else None)
    ev = load_anatuple(args.path, entry_stop=args.entry_stop, mc=mc)
    print(_summary(ev))
    if args.out:
        save_npz(ev, args.out)
        print(f"  saved -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
