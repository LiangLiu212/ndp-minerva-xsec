"""Truth-level observables as pure functions of a TruthTable (GeV, detector frame).

Every observable a channel or measurement manifest can name lives here, keyed by name in
`OBSERVABLES`. Adding an observable is adding a function and a registry entry. A manifest may
also give a numpy *expression* over the truth columns and the named observables (the same
namespace a `reweight` model's `weight_expr` sees), so a user can define an observable
without touching the code: `"sqrt(Q2)"`, `"E_nu - lep_E"`, `"lep_p*cos(lep_theta)"`.
"""
from __future__ import annotations

import keyword
import re

import numpy as np

from ..events import TruthTable, M_MU, M_P, M_N, INT_CODE

#: NuMI beam points ~3.34 degrees downward relative to the MINERvA detector z axis.
#: MAT's TruthFunctions::GetThetalepTrue rotates the lepton about x by this angle
#: (MinervaUnits::numi_beam_angle_rad). The channel decides the frame (`phase_space.frame`).
NUMI_BEAM_ANGLE_RAD = -0.05887

#: numpy names an expression may use (shared with the reco-observable expressions).
MATH_NAMESPACE = {"np": np, "where": np.where, "exp": np.exp, "log": np.log, "log10": np.log10, "sqrt": np.sqrt,
                  "abs": np.abs, "clip": np.clip, "minimum": np.minimum, "maximum": np.maximum, "pi": np.pi,
                  "sin": np.sin, "cos": np.cos, "tan": np.tan, "arccos": np.arccos, "arcsin": np.arcsin,
                  "arctan": np.arctan, "arctan2": np.arctan2, "deg2rad": np.deg2rad, "rad2deg": np.rad2deg,
                  "isfinite": np.isfinite, "nan": np.nan, "inf": np.inf}

_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def is_identifier(s: str) -> bool:
    return bool(_IDENT.match(s)) and not keyword.iskeyword(s)


def eval_expr(expr: str, ns: dict, n: int) -> np.ndarray:
    """Evaluate a manifest expression in a numpy-only namespace; broadcast scalars to n."""
    scope = dict(MATH_NAMESPACE)
    scope.update(ns)
    try:
        v = eval(expr, {"__builtins__": {}}, scope)  # noqa: S307 — analyst-authored, numpy namespace only
    except NameError as e:
        raise KeyError(f"expression {expr!r}: {e}; known names: {sorted(k for k in ns if not k.startswith('_'))}") from None
    return np.broadcast_to(np.asarray(v, float), (n,)).copy()


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


def lep_theta_deg(t: TruthTable, frame: str = "detector") -> np.ndarray:
    return np.rad2deg(lep_theta(t, frame))


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


def x_bj(t: TruthTable, **_) -> np.ndarray:
    """Bjorken x = Q2 / (2 M_N q0)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(q0(t) > 0, t["Q2"] / (2.0 * M_N * q0(t)), np.nan)


def y_inel(t: TruthTable, **_) -> np.ndarray:
    """Inelasticity y = q0 / E_nu."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(t["E_nu"] > 0, q0(t) / t["E_nu"], np.nan)


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


def E_had_fs(t: TruthTable, **_) -> np.ndarray:
    """Total final-state hadronic energy (every non-lepton particle, incl. neutrons) [GeV]."""
    if not t.has_fs:
        raise ValueError("E_had_fs needs final-state particles (fs_* columns)")
    apdg = np.abs(t["fs_pdg"])
    lep = np.isin(apdg, [11, 12, 13, 14, 15, 16]) | (apdg >= 1000000000)
    return t.fs_sum(np.where(lep, 0.0, t["fs_E"]))


def unit(t: TruthTable, **_) -> np.ndarray:
    """Constant 0.5: the dummy second axis of a one-dimensional measurement (edges [0, 1])."""
    return np.full(t.n, 0.5)


OBSERVABLES = {
    "lep_p": lep_p, "lep_pT": lep_pT, "lep_pz": lep_pz, "lep_theta": lep_theta, "lep_theta_deg": lep_theta_deg,
    "lep_E": lep_E, "E_nu": E_nu, "Q2": Q2, "W": W, "q0": q0, "q3": q3, "x_bj": x_bj, "y_inel": y_inel,
    "E_avail": E_avail, "E_had_fs": E_had_fs, "unit": unit,
}


def namespace(t: TruthTable, frame: str = "detector") -> dict:
    """Every truth column and every evaluable observable (in `frame`), plus the interaction codes."""
    ns = {k: t[k] for k in t.columns if k != "fs_offsets" and not k.startswith("fs_")}
    for name, f in OBSERVABLES.items():
        try:
            ns[name] = f(t, frame=frame)
        except Exception:  # e.g. E_avail on a sample without final-state particles
            pass
    ns.update(INT_CODE)          # QE, RES, DIS, COH, MEC as codes
    ns["M_MU"], ns["M_P"], ns["M_N"] = M_MU, M_P, M_N
    return ns


def evaluate(name_or_expr: str, t: TruthTable, **kw) -> np.ndarray:
    """A registered observable by name, or an expression over columns + observables."""
    s = str(name_or_expr).strip()
    if s in OBSERVABLES:
        return OBSERVABLES[s](t, **kw)
    if is_identifier(s):
        if s in t:
            return np.asarray(t[s], float)
        raise KeyError(f"unknown observable {s!r}; known: {sorted(OBSERVABLES)} or any truth column")
    return eval_expr(s, namespace(t, kw.get("frame", "detector")), t.n)
