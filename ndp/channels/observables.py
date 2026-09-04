"""Truth-level observables as pure functions of a TruthTable (GeV, detector frame).

Every observable a channel manifest can name lives here, keyed by name in
`OBSERVABLES`. Adding an observable is adding a function and a registry entry.
"""
from __future__ import annotations

import numpy as np

from ..events import TruthTable, M_MU, M_P, M_N

#: NuMI beam points ~3.34 degrees downward relative to the MINERvA detector z axis.
#: MAT's TruthFunctions::GetThetalepTrue rotates the lepton about x by this angle
#: (MinervaUnits::numi_beam_angle_rad). The platform default frame is *detector*
#: (matches the audited 2D run); pass frame="beam" to rotate.
NUMI_BEAM_ANGLE_RAD = -0.05887


def _lep_p3(t: TruthTable, frame: str = "detector"):
    px, py, pz = t["lep_px"], t["lep_py"], t["lep_pz"]
    if frame == "beam":
        a = NUMI_BEAM_ANGLE_RAD
        py, pz = py * np.cos(a) - pz * np.sin(a), py * np.sin(a) + pz * np.cos(a)
    elif frame != "detector":
        raise ValueError(f"unknown frame {frame!r}")
    return px, py, pz


def lep_p(t: TruthTable, frame: str = "detector") -> np.ndarray:
    px, py, pz = _lep_p3(t, frame)
    return np.sqrt(px * px + py * py + pz * pz)


def lep_pT(t: TruthTable, frame: str = "detector") -> np.ndarray:
    px, py, _ = _lep_p3(t, frame)
    return np.sqrt(px * px + py * py)


def lep_pz(t: TruthTable, frame: str = "detector") -> np.ndarray:
    return _lep_p3(t, frame)[2]


def lep_theta(t: TruthTable, frame: str = "detector") -> np.ndarray:
    p = lep_p(t, frame)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(p > 0, _lep_p3(t, frame)[2] / p, 1.0)
    return np.arccos(np.clip(c, -1.0, 1.0))


def lep_E(t: TruthTable, **_) -> np.ndarray:
    return t["lep_E"]


def E_nu(t: TruthTable, **_) -> np.ndarray:
    return t["E_nu"]


def Q2(t: TruthTable, **_) -> np.ndarray:
    return t["Q2"]


def W(t: TruthTable, **_) -> np.ndarray:
    return t["W"]


def q0(t: TruthTable, **_) -> np.ndarray:
    """Energy transfer E_nu - E_lep [GeV]."""
    return t["E_nu"] - t["lep_E"]


def q3(t: TruthTable, **_) -> np.ndarray:
    """Three-momentum transfer sqrt(Q2 + q0^2) [GeV]."""
    q = q0(t)
    return np.sqrt(np.maximum(t["Q2"] + q * q, 0.0))


def E_avail(t: TruthTable, **_) -> np.ndarray:
    """MINERvA 'available energy' (arXiv:2110.13372 Sec. 3), GeV.

    sum T_proton + sum T_pi+- + sum E_other over final-state particles, excluding
    neutrons, the primary lepton and all neutrinos; strange baryons contribute
    E - M_nucleon (they carry a nucleon mass that is not deposited).
    """
    if not t.has_fs:
        raise ValueError("E_avail needs final-state particles (fs_* columns)")
    pdg, E = t["fs_pdg"], t["fs_E"]
    apdg = np.abs(pdg)
    contrib = np.array(E, dtype=float)
    contrib = np.where(pdg == 2212, E - M_P, contrib)                    # proton KE
    contrib = np.where(apdg == 211, E - 0.13957039, contrib)             # charged pion KE
    strange_baryon = np.isin(apdg, [3122, 3222, 3212, 3112, 3322, 3312, 3334])
    contrib = np.where(strange_baryon, E - M_N, contrib)
    exclude = (pdg == 2112) | np.isin(apdg, [12, 14, 16]) | np.isin(apdg, [11, 13, 15])
    # nuclear remnants / pseudo-particles (GENIE codes > 1e9, 2000000xxx) deposit nothing
    exclude |= (apdg >= 1000000000) | ((apdg >= 2000000000) & (apdg < 3000000000))
    return t.fs_sum(np.where(exclude, 0.0, contrib))


OBSERVABLES = {
    "lep_p": lep_p, "lep_pT": lep_pT, "lep_pz": lep_pz, "lep_theta": lep_theta, "lep_E": lep_E,
    "E_nu": E_nu, "Q2": Q2, "W": W, "q0": q0, "q3": q3, "E_avail": E_avail,
}


def evaluate(name: str, t: TruthTable, **kw) -> np.ndarray:
    try:
        f = OBSERVABLES[name]
    except KeyError:
        raise KeyError(f"unknown observable {name!r}; known: {sorted(OBSERVABLES)}") from None
    return f(t, **kw)
