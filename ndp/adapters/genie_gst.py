"""GENIE `gst` summary ntuple -> TruthTable (GeV throughout; gst is already in GeV)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..events import TruthTable, Normalization

SCALARS = ["neu", "Ev", "fspl", "El", "pxl", "pyl", "pzl", "cc", "nc", "qel", "res", "dis", "coh", "mec",
           "Z", "A", "Q2", "W", "wght", "vtxx", "vtxy", "vtxz"]
FS = ["pdgf", "Ef", "pxf", "pyf", "pzf"]


def read_gst(path: str | Path, entry_stop: int | None = None, with_fs: bool = True,
             meta: dict | None = None) -> TruthTable:
    import uproot
    import awkward as ak
    tree = uproot.open(path)["gst"]
    a = tree.arrays(SCALARS, library="np", entry_stop=entry_stop)
    n = len(a["Ev"])
    int_type = np.zeros(n, dtype=np.int64)
    for code, key in ((1, "qel"), (2, "res"), (3, "dis"), (4, "coh"), (5, "mec")):
        int_type[a[key].astype(bool)] = code
    cols = {
        "nu_pdg": a["neu"], "E_nu": a["Ev"], "lep_pdg": a["fspl"],
        "lep_px": a["pxl"], "lep_py": a["pyl"], "lep_pz": a["pzl"], "lep_E": a["El"],
        "current": np.where(a["cc"].astype(bool), 1, 2), "int_type": int_type,
        "target_Z": a["Z"], "target_A": a["A"], "Q2": a["Q2"], "W": a["W"], "weight": a["wght"],
        # gst vertices are in metres for geometry runs and 0 for point targets; keep mm for symmetry
        "vtx_x": a["vtxx"] * 1e3, "vtx_y": a["vtxy"] * 1e3, "vtx_z": a["vtxz"] * 1e3,
    }
    if with_fs:
        parts = tree.arrays(FS, library="ak", entry_stop=entry_stop)
        counts = ak.num(parts["pdgf"], axis=1).to_numpy()
        cols["fs_offsets"] = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
        cols["fs_pdg"] = ak.flatten(parts["pdgf"]).to_numpy()
        cols["fs_E"] = ak.flatten(parts["Ef"]).to_numpy()
        cols["fs_px"] = ak.flatten(parts["pxf"]).to_numpy()
        cols["fs_py"] = ak.flatten(parts["pyf"]).to_numpy()
        cols["fs_pz"] = ak.flatten(parts["pzf"]).to_numpy()
    has_geom = bool(np.any(a["vtxz"] != 0.0))
    m = {"source": f"{path}:gst", "generator": "GENIE", "units": "GeV, GeV^2, mm",
         "n_generated": int(tree.num_entries), "norm": Normalization(kind="shape").to_dict(),
         "has_geometry": has_geom,
         "geometry_note": "point target: vertex phase-space cuts are not applied" if not has_geom else "vertices from a geometry run"}
    if meta:
        m.update(meta)
    return TruthTable(cols, m)
