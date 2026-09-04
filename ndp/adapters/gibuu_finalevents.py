"""GiBUU `FinalEvents.dat` adapter (neutrino mode, eventtype 5) -> TruthTable.

Column layout (GiBUU 2025, analysis/neutrinoAnalysis.f90):
  1 run  2 event  3 ID  4 charge  5 perweight  6-8 position [fm]  9-12 four-momentum (E,px,py,pz) [GeV]
  13 history  14 production_ID  15 E_nu [GeV]
Every particle row of an event carries the event's `perweight`; the struck nucleon is written with
weight 0 (and, with outputEvents_lepIn, so is the incoming lepton). The neutrino travels along +z.

production_ID (the first interaction): 1 QE; 2-31 baryon resonance (IdTable numbering); 32/33 pi
background off n/p; 34 DIS; 35 2p2h QE; 36 2p2h Delta; 37 two-pion background.
NDP mapping: 1->QE, 2-31->RES, 34->DIS, 35/36->MEC, 32/33/37->DIS (non-resonant pion production is
counted with DIS, the GENIE convention; the raw code is kept in `gibuu_production_id`).

Normalisation: perweights are in 1e-38 cm^2 per nucleon and sum, per run, to the total cross section,
so the per-event weight is perweight / n_runs and 1e-38 cm^2 per unit weight (the adapter cross-checks
against `neutrino_absorption_cross_section_ALL.dat` when it sits next to the file).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..events import TruthTable, Normalization, M_P, M_N

_MESON = {101: {1: 211, 0: 111, -1: -211}, 102: {0: 221}, 103: {1: 213, 0: 113, -1: -213}, 105: {0: 223},
          106: {0: 331}, 107: {0: 333}, 110: {1: 321, 0: 311}, 111: {-1: -321, 0: -311},
          112: {1: 411, 0: 421}, 113: {-1: -411, 0: -421}}
_LEPTON = {901: 11, 902: 13, 903: 15, 911: 12, 912: 14, 913: 16}
_BARYON = {32: {0: 3122}, 33: {1: 3222, 0: 3212, -1: 3112}, 2: {2: 2224, 1: 2214, 0: 2114, -1: 1114}}


def gibuu_pdg(gid: int, charge: int) -> int:
    aid, anti = abs(int(gid)), int(gid) < 0
    if aid == 1:
        pdg = 2212 if charge == 1 else 2112
    elif aid in _BARYON:
        pdg = _BARYON[aid].get(int(charge), 0)
    elif 3 <= aid <= 31:
        pdg = 2212 if charge == 1 else 2112     # undecayed N*: no PDG assignment kept; count as nucleon
    elif aid in _MESON:
        return _MESON[aid].get(int(charge), 0)  # mesons: charge already encodes the antiparticle
    elif aid in _LEPTON:
        pdg = _LEPTON[aid]
    elif aid == 999:
        return 22
    elif aid >= 1000:
        return 1000000000 + aid  # residual nucleus 1000+A' -> tagged as nuclear remnant
    else:
        return 0
    return -pdg if anti else pdg


def production_to_int_type(pid: int) -> int:
    if pid == 1:
        return 1
    if 2 <= pid <= 31:
        return 2
    if pid in (34, 32, 33, 37):
        return 3
    if pid in (35, 36):
        return 5
    return 0


def read_finalevents(path: str | Path, *, target_Z: int, target_A: int, n_runs: int | None = None,
                     nu_pdg_hint: int | None = None) -> TruthTable:
    path = Path(path)
    raw = np.loadtxt(path, comments="#", ndmin=2)
    if raw.shape[1] < 15:
        raise ValueError(f"{path}: expected 15 columns, got {raw.shape[1]}")
    run, evt, gid, chg = raw[:, 0].astype(int), raw[:, 1].astype(int), raw[:, 2].astype(int), raw[:, 3].astype(int)
    w, p4, prod, enu = raw[:, 4], raw[:, 8:12], raw[:, 13].astype(int), raw[:, 14]
    key = run * 10_000_000 + evt
    uniq, first = np.unique(key, return_index=True)
    order = np.argsort(key, kind="stable")
    key_s = key[order]
    bounds = np.searchsorted(key_s, uniq)
    bounds = np.append(bounds, len(key_s))
    n = len(uniq)
    cols = {k: np.zeros(n) for k in ("E_nu", "lep_px", "lep_py", "lep_pz", "lep_E", "Q2", "W", "weight")}
    icol = {k: np.zeros(n, dtype=np.int64) for k in ("nu_pdg", "lep_pdg", "current", "int_type", "gibuu_production_id")}
    fs_pdg, fs_E, fs_px, fs_py, fs_pz, offsets = [], [], [], [], [], [0]
    n_no_lepton = 0
    for i in range(n):
        idx = order[bounds[i]:bounds[i + 1]]
        g, c, ww, mom = gid[idx], chg[idx], w[idx], p4[idx]
        pdgs = np.array([gibuu_pdg(a, b) for a, b in zip(g, c)])
        is_lep = np.isin(np.abs(pdgs), [11, 13, 15, 12, 14, 16])
        out_lep = np.where(is_lep & (ww != 0))[0]
        struck = np.where((ww == 0) & (np.abs(pdgs) == 2212) | (ww == 0) & (np.abs(pdgs) == 2112))[0]
        e_nu = float(enu[idx][0])
        cols["E_nu"][i] = e_nu
        cols["weight"][i] = float(ww[ww != 0][0]) if np.any(ww != 0) else 0.0
        icol["gibuu_production_id"][i] = int(prod[idx][0])
        icol["int_type"][i] = production_to_int_type(int(prod[idx][0]))
        if out_lep.size:
            j = out_lep[0]
            lp = pdgs[j]
            cols["lep_E"][i], cols["lep_px"][i], cols["lep_py"][i], cols["lep_pz"][i] = mom[j]
            icol["lep_pdg"][i] = lp
            charged = abs(lp) in (11, 13, 15)
            icol["current"][i] = 1 if charged else 2
            flavour = abs(lp) + 1 if charged else abs(lp)
            icol["nu_pdg"][i] = int(np.sign(lp)) * flavour if charged else lp
            cols["Q2"][i] = max(2.0 * e_nu * (mom[j][0] - mom[j][3]) - (mom[j][0] ** 2 - np.sum(mom[j][1:] ** 2)), 0.0)
            if struck.size:
                pn = mom[struck[0]]
                tot = np.array([e_nu + pn[0] - mom[j][0], -mom[j][1] + pn[1], -mom[j][2] + pn[2], e_nu + pn[3] - mom[j][3]])
                cols["W"][i] = np.sqrt(max(tot[0] ** 2 - np.sum(tot[1:] ** 2), 0.0))
            else:
                q0 = e_nu - mom[j][0]
                cols["W"][i] = np.sqrt(max(M_N ** 2 + 2 * M_N * q0 - cols["Q2"][i], 0.0))
        else:
            n_no_lepton += 1
            icol["nu_pdg"][i] = nu_pdg_hint or 14
        keep = (ww != 0) & ~is_lep
        fs_pdg.extend(pdgs[keep]); fs_E.extend(mom[keep, 0]); fs_px.extend(mom[keep, 1]); fs_py.extend(mom[keep, 2]); fs_pz.extend(mom[keep, 3])
        offsets.append(len(fs_pdg))
    nr = int(n_runs or run.max())
    cols.update(icol)
    cols["target_Z"] = np.full(n, int(target_Z)); cols["target_A"] = np.full(n, int(target_A))
    cols["weight"] = cols["weight"] / nr
    cols.update({"fs_offsets": np.array(offsets, dtype=np.int64), "fs_pdg": np.array(fs_pdg, dtype=np.int64),
                 "fs_E": np.array(fs_E), "fs_px": np.array(fs_px), "fs_py": np.array(fs_py), "fs_pz": np.array(fs_pz)})
    meta = {"source": str(path), "generator": "GiBUU", "units": "GeV", "has_geometry": False, "n_runs": nr,
            "n_events_without_lepton_row": n_no_lepton,
            "norm": Normalization(kind="xsec_per_nucleon", xsec_per_unit_weight=1e-38,
                                  notes="perweight in 1e-38 cm^2/nucleon summing to sigma_tot per run; weights divided by n_runs").to_dict()}
    xs_file = path.parent / "neutrino_absorption_cross_section_ALL.dat"
    if xs_file.exists():
        try:
            last = [ln for ln in xs_file.read_text().splitlines() if ln.strip() and not ln.startswith("#")][-1]
            meta["gibuu_total_xsec_file_last_row"] = last.split()[:3]
            meta["sum_weights_1e-38cm2"] = float(cols["weight"].sum())
        except Exception:  # pragma: no cover
            pass
    return TruthTable(cols, meta)
