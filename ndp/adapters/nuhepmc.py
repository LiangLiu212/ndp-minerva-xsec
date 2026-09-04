"""NuHepMC (HepMC3) adapter — ACHILLES and any NuHepMC-compliant generator -> TruthTable.

Uses pyhepmc. Beam = status 4, target = status 20 (or 11), final state = status 1. The process
ID (event attribute `signal_process_id`, NuHepMC E.C.1) is mapped through the run's
`NuHepMC.ProcessInfo[<id>].Name` text (QE / MEC|2p2h / RES / DIS / COH) with the NuHepMC ranges
(2xx QE, 3xx MEC, 4xx RES, 5xx SIS, 6xx DIS, 7xx COH) as fallback. Cross sections: the
`NuHepMC.FluxAveragedTotalCrossSection` run attribute (G.C.4) when present, else the last event's
GenCrossSection; converted from the declared unit (pb / nb / ...) and target scale (PerAtom ->
divided by A) to cm^2 per nucleon.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..events import TruthTable, Normalization, M_N

_UNIT_TO_CM2 = {"pb": 1e-36, "nb": 1e-33, "fb": 1e-39, "mb": 1e-27, "cm2": 1.0}


def _attr(container, name, default=None):
    try:
        a = container.attributes[name]
    except Exception:
        return default
    try:
        return a.astype(str)
    except Exception:
        return str(a)


def name_to_int_type(name: str, pid: int) -> int:
    n = (name or "").lower()
    for key, code in (("coh", 4), ("mec", 5), ("2p2h", 5), ("qe", 1), ("quasi", 1), ("res", 2), ("dis", 3)):
        if key in n:
            return code
    return {2: 1, 3: 5, 4: 2, 5: 3, 6: 3, 7: 4}.get(pid // 100, 0)


def read_nuhepmc(path: str | Path, entry_stop: int | None = None) -> TruthTable:
    import pyhepmc
    rows = {k: [] for k in ("nu_pdg", "E_nu", "lep_pdg", "lep_px", "lep_py", "lep_pz", "lep_E", "current", "int_type",
                            "target_Z", "target_A", "Q2", "W", "weight", "process_id")}
    fs = {k: [] for k in ("pdg", "E", "px", "py", "pz")}
    offsets = [0]
    run_attrs, proc_names, last_xs, n_read = {}, {}, None, 0
    with pyhepmc.open(str(path)) as f:
        for evt in f:
            if n_read == 0 and evt.run_info is not None:
                for k in evt.run_info.attributes:
                    run_attrs[k] = _attr(evt.run_info, k, "")
                for k, v in run_attrs.items():
                    if k.startswith("NuHepMC.ProcessInfo[") and k.endswith("].Name"):
                        proc_names[int(k[len("NuHepMC.ProcessInfo["):-len("].Name")])] = v
            parts = list(evt.particles)
            beam = [p for p in parts if p.status == 4]
            tgt = [p for p in parts if p.status in (20, 11)]
            final = [p for p in parts if p.status == 1]
            nu = next((p for p in beam if abs(p.pid) in (12, 14, 16, 11)), beam[0] if beam else None)
            if nu is None:
                continue
            nupid = nu.pid
            lep = next((p for p in final if abs(p.pid) in (11, 13, 15) and np.sign(p.pid) == np.sign(nupid)), None)
            current = 1
            if lep is None:
                lep = next((p for p in final if p.pid == nupid), None)
                current = 2
            pid_attr = _attr(evt, "signal_process_id", None) or _attr(evt, "ProcID", "0")
            pid = int(float(pid_attr)) if pid_attr else 0
            tpdg = tgt[0].pid if tgt else 0
            A = (tpdg // 10) % 1000 if tpdg > 1000000000 else (1 if tpdg in (2212, 2112) else 0)
            Z = (tpdg // 10000) % 1000 if tpdg > 1000000000 else (1 if tpdg == 2212 else 0)
            k = nu.momentum
            rows["nu_pdg"].append(nupid); rows["E_nu"].append(k.e)
            if lep is not None:
                kp = lep.momentum
                q0, qx, qy, qz = k.e - kp.e, k.px - kp.px, k.py - kp.py, k.pz - kp.pz
                Q2 = max(qx * qx + qy * qy + qz * qz - q0 * q0, 0.0)
                rows["lep_pdg"].append(lep.pid); rows["lep_px"].append(kp.px); rows["lep_py"].append(kp.py); rows["lep_pz"].append(kp.pz); rows["lep_E"].append(kp.e)
                rows["Q2"].append(Q2); rows["W"].append(np.sqrt(max(M_N ** 2 + 2 * M_N * q0 - Q2, 0.0)))
            else:
                for kk in ("lep_pdg", "lep_px", "lep_py", "lep_pz", "lep_E", "Q2", "W"):
                    rows[kk].append(0)
            rows["current"].append(current); rows["int_type"].append(name_to_int_type(proc_names.get(pid, ""), pid))
            rows["target_Z"].append(Z); rows["target_A"].append(A); rows["process_id"].append(pid)
            rows["weight"].append(float(evt.weights[0]) if len(evt.weights) else 1.0)
            for p in final:
                if p is lep:
                    continue
                fs["pdg"].append(p.pid); fs["E"].append(p.momentum.e); fs["px"].append(p.momentum.px); fs["py"].append(p.momentum.py); fs["pz"].append(p.momentum.pz)
            offsets.append(len(fs["pdg"]))
            try:
                xs = evt.cross_section
                if xs is not None:
                    last_xs = float(xs.xsec(0))
            except Exception:
                pass
            n_read += 1
            if entry_stop and n_read >= entry_stop:
                break
    cols = {k: np.asarray(v) for k, v in rows.items()}
    cols.update({"fs_offsets": np.asarray(offsets, dtype=np.int64), "fs_pdg": np.asarray(fs["pdg"], dtype=np.int64),
                 "fs_E": np.asarray(fs["E"]), "fs_px": np.asarray(fs["px"]), "fs_py": np.asarray(fs["py"]), "fs_pz": np.asarray(fs["pz"])})
    unit = run_attrs.get("NuHepMC.Units.CrossSection.Unit", "pb")
    scale = run_attrs.get("NuHepMC.Units.CrossSection.TargetScale", "PerAtom")
    xs_total = run_attrs.get("NuHepMC.FluxAveragedTotalCrossSection")
    xs_val = float(xs_total) if xs_total else last_xs
    A_mode = int(np.median(cols["target_A"])) if n_read else 1
    norm = Normalization(kind="shape")
    if xs_val is not None and n_read:
        per_nucleon = xs_val * _UNIT_TO_CM2.get(unit, 1e-36) / (A_mode if scale == "PerAtom" and A_mode else 1.0)
        norm = Normalization(kind="xsec_per_nucleon", xsec_per_unit_weight=per_nucleon / float(cols["weight"].sum()),
                             notes=f"{'G.C.4 flux-averaged total' if xs_total else 'last-event GenCrossSection'} {xs_val:g} {unit} {scale} -> cm^2/nucleon")
    meta = {"source": str(path), "generator": run_attrs.get("NuHepMC.Generator", "NuHepMC generator"), "units": "GeV",
            "has_geometry": False, "n_generated": n_read, "process_names": proc_names, "norm": norm.to_dict(),
            "run_attributes": {k: v for k, v in run_attrs.items() if not k.startswith("NuHepMC.ProcessInfo")}}
    return TruthTable(cols, meta)
