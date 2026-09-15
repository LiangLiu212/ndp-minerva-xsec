"""Truth-level observables as pure functions of a TruthTable (GeV, detector frame).

Every observable a channel or measurement manifest can name lives here, keyed by name in
`OBSERVABLES`. Adding an observable is adding a function and a registry entry. A manifest may
also give a numpy *expression* over the truth columns and the named observables (the same
namespace a `reweight` model's `weight_expr` sees), so a user can define an observable
without touching the code: `"sqrt(Q2)"`, `"E_nu - lep_E"`, `"lep_p*cos(lep_theta)"`.
"""
from __future__ import annotations

import json
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


class LazyNamespace(dict):
    """A dict of columns/constants whose missing names are computed on first use by `resolvers`
    ({name: callable}). Lets `eval` see every observable without evaluating every observable."""

    def __init__(self, base: dict, resolvers: dict):
        super().__init__(base)
        self._resolvers = resolvers

    def __missing__(self, key):
        f = self._resolvers.get(key)
        if f is None:
            raise KeyError(key)
        v = f()
        self[key] = v
        return v

    def known_names(self) -> list:
        return sorted(set(self) | set(self._resolvers))


def eval_expr(expr: str, ns: dict, n: int) -> np.ndarray:
    """Evaluate a manifest expression in a numpy-only namespace; broadcast scalars to n."""
    if isinstance(ns, LazyNamespace):
        scope = LazyNamespace(dict(MATH_NAMESPACE, **{k: ns[k] for k in dict.keys(ns)}), ns._resolvers)
        known = ns.known_names()
    else:
        scope = dict(MATH_NAMESPACE)
        scope.update(ns)
        known = sorted(k for k in ns if not k.startswith("_"))
    try:
        v = eval(expr, {"__builtins__": {}}, scope)  # noqa: S307 — analyst-authored, numpy namespace only
    except (NameError, KeyError) as e:
        raise KeyError(f"expression {expr!r}: {e}; known names: {known}") from None
    return np.broadcast_to(np.asarray(v, float), (n,)).copy()


FRAMES = ("detector", "beam")


def native_frame(t: TruthTable) -> str:
    """The frame a table's momenta are stored in.

    `meta["frame"]` when the producer recorded it (the MINERvA adapter writes "detector"); otherwise a
    generator sample without detector geometry fires its neutrino along +z, which *is* the beam axis,
    so it is "beam"; anything else (toy tables, legacy caches) is taken as "detector", the historical
    behaviour. Before 2026-09-15 every table was rotated as if it were detector-frame, which tilted
    generator samples by the 3.4 degree NuMI angle in beam-frame channels (docs/decisions.md).
    """
    m = getattr(t, "meta", None) or {}
    f = m.get("frame")
    if f in FRAMES:
        return f
    if m.get("generator") and not m.get("has_geometry", False):
        return "beam"
    return "detector"


def frame_rotation_angle(native: str, requested: str) -> float:
    """Rotation about x (radians) that takes `native`-frame momenta into the `requested` frame."""
    for f in (native, requested):
        if f not in FRAMES:
            raise ValueError(f"unknown frame {f!r} (known: {FRAMES})")
    if native == requested:
        return 0.0
    return NUMI_BEAM_ANGLE_RAD if (native, requested) == ("detector", "beam") else -NUMI_BEAM_ANGLE_RAD


def rotate_about_x(px, py, pz, angle: float):
    if angle == 0.0:
        return px, py, pz
    c, s = np.cos(angle), np.sin(angle)
    return px, py * c - pz * s, py * s + pz * c


def _lep_p3(t: TruthTable, frame: str = "detector"):
    px, py, pz = t["lep_px"], t["lep_py"], t["lep_pz"]
    return rotate_about_x(px, py, pz, frame_rotation_angle(native_frame(t), frame))


def lep_p(t: TruthTable, frame: str = "detector", **_) -> np.ndarray:
    px, py, pz = _lep_p3(t, frame)
    return np.sqrt(px * px + py * py + pz * pz)


def lep_pT(t: TruthTable, frame: str = "detector", **_) -> np.ndarray:
    px, py, _ = _lep_p3(t, frame)
    return np.sqrt(px * px + py * py)


def lep_pz(t: TruthTable, frame: str = "detector", **_) -> np.ndarray:
    return _lep_p3(t, frame)[2]


def lep_theta(t: TruthTable, frame: str = "detector", **_) -> np.ndarray:
    p = lep_p(t, frame)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(p > 0, _lep_p3(t, frame)[2] / p, 1.0)
    return np.arccos(np.clip(c, -1.0, 1.0))


def lep_theta_deg(t: TruthTable, frame: str = "detector", **_) -> np.ndarray:
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
        if "E_avail" in t:                      # a skim stores the derived value
            return np.asarray(t["E_avail"], float)
        raise ValueError("E_avail needs final-state particles (fs_* columns) or a derived E_avail column")
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
        if "E_had_fs" in t:
            return np.asarray(t["E_had_fs"], float)
        raise ValueError("E_had_fs needs final-state particles (fs_* columns) or a derived E_had_fs column")
    apdg = np.abs(t["fs_pdg"])
    lep = np.isin(apdg, [11, 12, 13, 14, 15, 16]) | (apdg >= 1000000000)
    return t.fs_sum(np.where(lep, 0.0, t["fs_E"]))


def unit(t: TruthTable, **_) -> np.ndarray:
    """Constant 0.5: the dummy second axis of a one-dimensional measurement (edges [0, 1])."""
    return np.full(t.n, 0.5)


# ---- leading proton (the channel's signal block says which protons count) -----------------
def _leading(t: TruthTable, frame: str, signal: dict | None):
    from .signal import leading_proton  # local import: signal.py imports this module
    window = (signal or {}).get("proton") if isinstance(signal, dict) else None
    return leading_proton(t, frame, window)


def proton_p(t: TruthTable, frame: str = "detector", signal: dict | None = None, **_) -> np.ndarray:
    """Momentum [GeV/c] of the leading proton inside the signal block's proton window (NaN if none)."""
    return _leading(t, frame, signal)["p"]


def proton_theta(t: TruthTable, frame: str = "detector", signal: dict | None = None, **_) -> np.ndarray:
    return _leading(t, frame, signal)["theta"]


def proton_theta_deg(t: TruthTable, frame: str = "detector", signal: dict | None = None, **_) -> np.ndarray:
    return np.rad2deg(_leading(t, frame, signal)["theta"])


def proton_pT(t: TruthTable, frame: str = "detector", signal: dict | None = None, **_) -> np.ndarray:
    return _leading(t, frame, signal)["pT"]


def n_protons_in_window(t: TruthTable, frame: str = "detector", signal: dict | None = None, **_) -> np.ndarray:
    return _leading(t, frame, signal)["n_in_window"].astype(float)


# ---- transverse kinematic imbalance (muon + leading proton) --------------------------------
# Definitions follow Lu et al., PRC 94 (2016) 015503 and Furmanski & Sobczyk, PRC 95 (2017) 065501,
# as used by MINERvA (arXiv:1805.05486 Eqs. 1-6; arXiv:2503.15047 Fig. 1 and Eqs. pl/pn). With
# z = neutrino direction (the frame's z axis after the beam rotation), pT the transverse momenta:
#     dpT_vec  = pT_mu + pT_p                       missing transverse momentum
#     dpT      = |dpT_vec|
#     dalphaT  = arccos( -pT_mu . dpT_vec / (|pT_mu| dpT) )          angle between -pT_mu and dpT_vec
#     dphiT    = arccos( -pT_mu . pT_p  / (|pT_mu| |pT_p|) )         deviation from back-to-back ("coplanarity")
#     dpTx     = (z x pT_mu_hat) . dpT_vec         component perpendicular to the muon transverse direction
#     dpTy     = -pT_mu_hat . dpT_vec              component along -pT_mu (negative: the proton carries less pT)
#     dpL      = R/2 - (m_A'^2 + dpT^2)/(2R),  R = m_A + pL_mu + pL_p - E_mu - E_p
#     pn       = sqrt(dpT^2 + dpL^2)               inferred initial-state neutron momentum
# m_A (target nucleus) and m_A' (residual nucleus, = m_A - m_n + b) come from the channel manifest's
# `observable_params.tki` block, never from code. MAT's legacy MnvRecoShifter::Calc_tki_vars uses the
# opposite sign for both dpTx and dpTy (its x axis is pT_mu x z, its y axis along +pT_mu); the
# released dpTy binning (tail to -6.5 GeV/c) is the Lu convention used here.
_TKI_CACHE = {}


def _tki(t: TruthTable, frame: str, signal: dict | None, params: dict | None) -> dict:
    from .signal import leading_proton, rotate_to_frame
    window = (signal or {}).get("proton") if isinstance(signal, dict) else None
    key = (id(t), frame, json.dumps(window, sort_keys=True))
    hit = _TKI_CACHE.get(key)
    if hit is not None and hit[0] is t:
        return hit[1]
    lp = leading_proton(t, frame, window)
    mx, my, mz = _lep_p3(t, frame)
    mE = t["lep_E"]
    px, py, pz, pE = lp["px"], lp["py"], lp["pz"], lp["E"]
    mT = np.sqrt(mx * mx + my * my)
    pT = np.sqrt(px * px + py * py)
    dx, dy = mx + px, my + py
    dpt = np.sqrt(dx * dx + dy * dy)
    with np.errstate(invalid="ignore", divide="ignore"):
        ux, uy = mx / mT, my / mT                                  # pT_mu_hat
        cos_alpha = -(ux * dx + uy * dy) / dpt
        cos_phi = -(ux * px + uy * py) / pT
        dalpha = np.where(dpt > 0, np.arccos(np.clip(cos_alpha, -1.0, 1.0)), np.nan)
        dphi = np.where(pT > 0, np.arccos(np.clip(cos_phi, -1.0, 1.0)), np.nan)
    dptx = -uy * dx + ux * dy                                      # (z x pT_mu_hat) . dpT_vec
    dpty = -(ux * dx + uy * dy)                                    # -pT_mu_hat . dpT_vec
    out = {"dpT": dpt, "dpTx": dptx, "dpTy": dpty, "dalphaT": dalpha, "dphiT": dphi,
           "_mz": mz, "_mE": mE, "_pz": pz, "_pE": pE}
    _TKI_CACHE.clear()
    _TKI_CACHE[key] = (t, out)
    return out


def _tki_masses(params: dict | None) -> tuple[float, float]:
    tki = (params or {}).get("tki") if isinstance(params, dict) else None
    if not tki or "m_A_gev" not in tki:
        raise KeyError("dpL / pn need the channel's `observable_params.tki` block (m_A_gev and m_Aprime_gev, "
                       "or m_A_gev + m_n_gev + excitation_b_gev); the platform carries no nuclear-mass defaults")
    m_A = float(tki["m_A_gev"])
    if "m_Aprime_gev" in tki:
        m_Ap = float(tki["m_Aprime_gev"])
    else:
        m_Ap = m_A - float(tki.get("m_n_gev", M_N)) + float(tki["excitation_b_gev"])
    return m_A, m_Ap


def dpT(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    """|pT_mu + pT_p| [GeV/c] with the leading proton of the signal block's window (NaN if none)."""
    return _tki(t, frame, signal, params)["dpT"]


def dpTx(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    return _tki(t, frame, signal, params)["dpTx"]


def dpTy(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    return _tki(t, frame, signal, params)["dpTy"]


def dalphaT(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    """Transverse boosting angle [rad]; NaN when dpT == 0."""
    return _tki(t, frame, signal, params)["dalphaT"]


def dalphaT_deg(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    return np.rad2deg(_tki(t, frame, signal, params)["dalphaT"])


def dphiT(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    """Coplanarity angle [rad]: deviation of the proton from back-to-back with the muon in the transverse plane."""
    return _tki(t, frame, signal, params)["dphiT"]


def dphiT_deg(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    return np.rad2deg(_tki(t, frame, signal, params)["dphiT"])


def dpL(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    """Longitudinal momentum imbalance [GeV/c] under the one-nucleon-knockout hypothesis (needs observable_params.tki)."""
    m_A, m_Ap = _tki_masses(params)
    k = _tki(t, frame, signal, params)
    R = m_A + k["_mz"] + k["_pz"] - k["_mE"] - k["_pE"]
    with np.errstate(invalid="ignore", divide="ignore"):
        return 0.5 * R - (m_Ap * m_Ap + k["dpT"] ** 2) / (2.0 * R)


def pn(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None, **_) -> np.ndarray:
    """Inferred initial-state neutron momentum sqrt(dpT^2 + dpL^2) [GeV/c]."""
    k = _tki(t, frame, signal, params)
    return np.sqrt(k["dpT"] ** 2 + dpL(t, frame, signal, params) ** 2)


OBSERVABLES = {
    "lep_p": lep_p, "lep_pT": lep_pT, "lep_pz": lep_pz, "lep_theta": lep_theta, "lep_theta_deg": lep_theta_deg,
    "lep_E": lep_E, "E_nu": E_nu, "Q2": Q2, "W": W, "q0": q0, "q3": q3, "x_bj": x_bj, "y_inel": y_inel,
    "E_avail": E_avail, "E_had_fs": E_had_fs, "unit": unit,
    "proton_p": proton_p, "proton_theta": proton_theta, "proton_theta_deg": proton_theta_deg, "proton_pT": proton_pT,
    "n_protons_in_window": n_protons_in_window,
    "dpT": dpT, "dpTx": dpTx, "dpTy": dpTy, "dalphaT": dalphaT, "dalphaT_deg": dalphaT_deg,
    "dphiT": dphiT, "dphiT_deg": dphiT_deg, "dpL": dpL, "pn": pn,
}


def namespace(t: TruthTable, frame: str = "detector", signal: dict | None = None, params: dict | None = None) -> dict:
    """Every truth column and every observable (in `frame`), plus the interaction codes. Observables
    are evaluated lazily, on first use by an expression (an unevaluable one, e.g. E_avail on a sample
    without final-state particles or dpL without tki masses, only fails if the expression names it)."""
    ns = {k: t[k] for k in t.columns if k != "fs_offsets" and not k.startswith("fs_")}
    ns.update(INT_CODE)          # QE, RES, DIS, COH, MEC as codes
    ns["M_MU"], ns["M_P"], ns["M_N"] = M_MU, M_P, M_N
    resolvers = {name: (lambda f=f: f(t, frame=frame, signal=signal, params=params)) for name, f in OBSERVABLES.items()}
    return LazyNamespace(ns, resolvers)


def evaluate(name_or_expr: str, t: TruthTable, **kw) -> np.ndarray:
    """A registered observable by name, or an expression over columns + observables."""
    s = str(name_or_expr).strip()
    if s in OBSERVABLES:
        return OBSERVABLES[s](t, **kw)
    if is_identifier(s):
        if s in t:
            return np.asarray(t[s], float)
        raise KeyError(f"unknown observable {s!r}; known: {sorted(OBSERVABLES)} or any truth column")
    return eval_expr(s, namespace(t, kw.get("frame", "detector"), kw.get("signal"), kw.get("params")), t.n)
