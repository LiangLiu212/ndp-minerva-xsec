"""NuWro `treeout` adapter (PyROOT + NuWro's event1 dictionary) -> TruthTable.

Needs `import ROOT` and `$NUWRO/bin/event1.so` (the pixi environment has both). Per event:
`in[0]` neutrino, `in[1]` struck nucleon (absent for coherent), `out` primary vertex products,
`post` final state after the cascade; `flag.{qel,res,dis,coh,mec}`, `flag.cc`; `weight` = the
total cross section in cm^2 (NuWro stores sigma_tot on every event when writing the file; events
are unweighted). NuWro quotes cross sections per nucleon — recorded as an assumption in docs.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from ..events import TruthTable, Normalization

MEV = 1e-3


def read_nuwro(path: str | Path, entry_stop: int | None = None, event1_so: str | None = None) -> TruthTable:
    import ROOT
    so = event1_so or os.path.join(os.environ.get("NUWRO", ""), "bin", "event1.so")
    if ROOT.gSystem.Load(so) < 0:
        raise RuntimeError(f"cannot load NuWro dictionary {so}")
    f = ROOT.TFile.Open(str(path))
    t = f.Get("treeout")
    rows = {k: [] for k in ("nu_pdg", "E_nu", "lep_pdg", "lep_px", "lep_py", "lep_pz", "lep_E", "current", "int_type",
                            "target_Z", "target_A", "Q2", "W", "weight", "nuwro_dyn")}
    fs = {k: [] for k in ("pdg", "E", "px", "py", "pz")}
    offsets, xsecs = [0], []
    for i, ev in enumerate(t):
        e = ev.e
        nu = e.in_[0] if hasattr(e, "in_") else getattr(e, "in")[0]
        lep = e.out[0]
        rows["nu_pdg"].append(int(nu.pdg)); rows["E_nu"].append(nu.t * MEV)
        rows["lep_pdg"].append(int(lep.pdg)); rows["lep_px"].append(lep.x * MEV); rows["lep_py"].append(lep.y * MEV)
        rows["lep_pz"].append(lep.z * MEV); rows["lep_E"].append(lep.t * MEV)
        rows["current"].append(2 if e.flag.nc else 1)
        it = 1 if e.flag.qel else 2 if e.flag.res else 3 if e.flag.dis else 4 if e.flag.coh else 5 if e.flag.mec else 0
        rows["int_type"].append(it); rows["nuwro_dyn"].append(int(e.dyn))
        rows["target_Z"].append(int(e.par.nucleus_p)); rows["target_A"].append(int(e.par.nucleus_p + e.par.nucleus_n))
        rows["Q2"].append(e.q2() * MEV * MEV if hasattr(e, "q2") else 0.0)
        rows["W"].append(e.W() * MEV if hasattr(e, "W") else 0.0)
        rows["weight"].append(1.0); xsecs.append(float(e.weight))
        for p in e.post:
            fs["pdg"].append(int(p.pdg)); fs["E"].append(p.t * MEV); fs["px"].append(p.x * MEV); fs["py"].append(p.y * MEV); fs["pz"].append(p.z * MEV)
        offsets.append(len(fs["pdg"]))
        if entry_stop and i + 1 >= entry_stop:
            break
    f.Close()
    cols = {k: np.asarray(v) for k, v in rows.items()}
    cols.update({"fs_offsets": np.asarray(offsets, dtype=np.int64), "fs_pdg": np.asarray(fs["pdg"], dtype=np.int64),
                 "fs_E": np.asarray(fs["E"]), "fs_px": np.asarray(fs["px"]), "fs_py": np.asarray(fs["py"]), "fs_pz": np.asarray(fs["pz"])})
    n = len(cols["E_nu"])
    sigma = float(np.mean(xsecs)) if xsecs else None
    meta = {"source": str(path), "generator": "NuWro", "units": "GeV", "has_geometry": False, "n_generated": n,
            "sigma_tot_cm2_per_nucleon": sigma,
            "norm": Normalization(kind="xsec_per_nucleon", xsec_per_unit_weight=sigma / n,
                                  notes="NuWro event.weight = sigma_tot [cm^2], assumed per nucleon (see open questions)").to_dict() if sigma else Normalization().to_dict()}
    return TruthTable(cols, meta)
