"""Reconstruction-level observables: functions of the cached reco columns of an experiment.

The reco side of a measurement is evaluated on a *reco table* — the dict of per-candidate
arrays the experiment adapter caches (`ndp data cache`), all in GeV / GeV^2 / mm. The
registry below names the derived quantities an analyst is likely to want; anything else
can be written as an expression over the cached columns and these names, e.g.
`"reco_E_mu + reco_recoil_E_calo"` or `"rad2deg(reco_theta)"`.

Cached columns (MINERvA MasterAnaDev adapter; see `adapters/minerva_anatuple.py::RECO_CACHE_COLUMNS`):
    reco_p reco_theta reco_pT reco_pz reco_E_mu reco_thetaX reco_thetaY reco_minos_p reco_minos_qp
    reco_recoil_E (MasterAnaDev_recoil_E)             reco_recoil_E_passive (…_recoil_passivecorrected)
    reco_recoil_E_calo (blob_recoil_E)                reco_recoil_E_polyline (recoil_E_polylinecorrected)
    reco_anatool_E_nu (MasterAnaDev_E)  reco_anatool_W  reco_anatool_x  reco_anatool_y
    reco_visible_E  reco_n_prongs  reco_n_hadron_tracks  reco_vtx_x reco_vtx_y reco_vtx_z
"""
from __future__ import annotations

import numpy as np

from ..events import M_MU
from .observables import eval_expr, is_identifier


def _col(r: dict, name: str) -> np.ndarray:
    try:
        return np.asarray(r[name], float)
    except KeyError:
        raise KeyError(f"reco column {name!r} is not in the cache (columns: {sorted(r)}); "
                       "rebuild it with `python -m ndp data cache`") from None


def reco_p(r, **_): return _col(r, "reco_p")
def reco_theta(r, **_): return _col(r, "reco_theta")
def reco_theta_deg(r, **_): return np.rad2deg(_col(r, "reco_theta"))
def reco_pT(r, **_): return _col(r, "reco_pT")
def reco_pz(r, **_): return _col(r, "reco_pz")
def reco_E_mu(r, **_): return _col(r, "reco_E_mu")
def reco_recoil_E(r, **_): return _col(r, "reco_recoil_E")
def reco_recoil_E_calo(r, **_): return _col(r, "reco_recoil_E_calo")
def reco_recoil_E_passive(r, **_): return _col(r, "reco_recoil_E_passive")
def reco_anatool_E_nu(r, **_): return _col(r, "reco_anatool_E_nu")
def reco_anatool_W(r, **_): return _col(r, "reco_anatool_W")
def reco_visible_E(r, **_): return _col(r, "reco_visible_E")
def reco_n_prongs(r, **_): return _col(r, "reco_n_prongs")
def reco_n_hadron_tracks(r, **_): return _col(r, "reco_n_hadron_tracks")


def reco_E_nu_calo(r, **_):
    """Calorimetric neutrino energy: muon energy + MasterAnaDev_recoil_E (the CC-inclusive recoil)."""
    return reco_E_mu(r) + reco_recoil_E(r)


def reco_q0_calo(r, **_):
    return reco_recoil_E(r)


def reco_Q2_calo(r, **_):
    """Q^2 = 2 E_nu (E_mu - p_mu cos theta) - m_mu^2 with the calorimetric E_nu [GeV^2]."""
    return 2.0 * reco_E_nu_calo(r) * (reco_E_mu(r) - reco_p(r) * np.cos(reco_theta(r))) - M_MU ** 2


def reco_q3_calo(r, **_):
    q0 = reco_q0_calo(r)
    return np.sqrt(np.maximum(reco_Q2_calo(r) + q0 * q0, 0.0))


def reco_Q2_anatool(r, **_):
    """Q^2 from the analysis tool's own E_nu (MasterAnaDev_E) and the muon [GeV^2]."""
    return 2.0 * reco_anatool_E_nu(r) * (reco_E_mu(r) - reco_p(r) * np.cos(reco_theta(r))) - M_MU ** 2


def unit(r, **_):
    n = len(next(iter(v for k, v in r.items() if k != "__meta__")))
    return np.full(n, 0.5)


# ---- leading proton candidate and transverse kinematic imbalance (cache v3) ---------------
# The cached MasterAnaDev proton candidate (momentum by dE/dx range) is the tool's primary candidate,
# which is the highest-momentum one on the open-data MC. Its momentum components and the muon's
# (MasterAnaDev_leptonE) are detector-frame; both are rotated into the NuMI beam frame here, the
# frame the tuple's own angle branches (muon_thetaX/Y, MasterAnaDev_proton_theta) use, so the reco
# TKI variables have exactly the truth-level definitions of `observables.py` (Lu et al. conventions).
_BEAM = -0.05887  # NUMI_BEAM_ANGLE_RAD (observables.py)


def _rot(px, py, pz):
    a = _BEAM
    return px, py * np.cos(a) - pz * np.sin(a), py * np.sin(a) + pz * np.cos(a)


def _tki(r: dict, params: dict | None = None) -> dict:
    mx, my, mz = _rot(_col(r, "reco_mu_px"), _col(r, "reco_mu_py"), _col(r, "reco_mu_pz"))
    px, py, pz = _rot(_col(r, "reco_proton_px"), _col(r, "reco_proton_py"), _col(r, "reco_proton_pz"))
    mT, pT = np.hypot(mx, my), np.hypot(px, py)
    dx, dy = mx + px, my + py
    dpt = np.hypot(dx, dy)
    with np.errstate(invalid="ignore", divide="ignore"):
        ux, uy = mx / mT, my / mT
        dalpha = np.where(dpt > 0, np.arccos(np.clip(-(ux * dx + uy * dy) / dpt, -1, 1)), np.nan)
        dphi = np.where(pT > 0, np.arccos(np.clip(-(ux * px + uy * py) / pT, -1, 1)), np.nan)
    out = {"dpT": dpt, "dpTx": -uy * dx + ux * dy, "dpTy": -(ux * dx + uy * dy), "dalphaT": dalpha, "dphiT": dphi,
           "proton_pT": pT, "proton_p": np.sqrt(px * px + py * py + pz * pz),
           "proton_theta": np.arccos(np.clip(np.where(pT + np.abs(pz) > 0, pz / np.sqrt(px * px + py * py + pz * pz), 1.0), -1, 1))}
    tki = (params or {}).get("tki") if isinstance(params, dict) else None
    if tki and "m_A_gev" in tki:
        m_A = float(tki["m_A_gev"])
        m_Ap = float(tki["m_Aprime_gev"]) if "m_Aprime_gev" in tki else m_A - float(tki.get("m_n_gev", 0.93956542)) + float(tki["excitation_b_gev"])
        R = m_A + mz + pz - _col(r, "reco_E_mu") - _col(r, "reco_proton_E")
        with np.errstate(invalid="ignore", divide="ignore"):
            out["dpL"] = 0.5 * R - (m_Ap * m_Ap + dpt * dpt) / (2.0 * R)
        out["pn"] = np.sqrt(dpt * dpt + out["dpL"] ** 2)
    return out


def _need_masses(k: dict, name: str):
    if name not in k:
        raise KeyError(f"{name} needs the channel's observable_params.tki block (m_A_gev, ...); pass params=channel.observable_params")
    return k[name]


def reco_proton_p(r, **_): return _col(r, "reco_proton_p")
def reco_proton_T(r, **_): return _col(r, "reco_proton_T")
def reco_proton_theta(r, **_): return _col(r, "reco_proton_theta_beam")          # tuple's beam-frame angle
def reco_proton_theta_deg(r, **_): return np.rad2deg(_col(r, "reco_proton_theta_beam"))
def reco_proton_pT(r, params=None, **_): return _tki(r, params)["proton_pT"]
def reco_proton_score1(r, **_): return _col(r, "reco_proton_score1")
def reco_n_michel(r, **_): return _col(r, "reco_n_michel")
def reco_n_iso_blobs(r, **_): return _col(r, "reco_n_iso_blobs")
def reco_dpT(r, params=None, **_): return _tki(r, params)["dpT"]
def reco_dpTx(r, params=None, **_): return _tki(r, params)["dpTx"]
def reco_dpTy(r, params=None, **_): return _tki(r, params)["dpTy"]
def reco_dalphaT(r, params=None, **_): return _tki(r, params)["dalphaT"]
def reco_dalphaT_deg(r, params=None, **_): return np.rad2deg(_tki(r, params)["dalphaT"])
def reco_dphiT(r, params=None, **_): return _tki(r, params)["dphiT"]
def reco_dphiT_deg(r, params=None, **_): return np.rad2deg(_tki(r, params)["dphiT"])
def reco_dpL(r, params=None, **_): return _need_masses(_tki(r, params), "dpL")
def reco_pn(r, params=None, **_): return _need_masses(_tki(r, params), "pn")


RECO_OBSERVABLES = {
    "reco_p": reco_p, "reco_theta": reco_theta, "reco_theta_deg": reco_theta_deg, "reco_pT": reco_pT, "reco_pz": reco_pz,
    "reco_E_mu": reco_E_mu, "reco_recoil_E": reco_recoil_E, "reco_recoil_E_calo": reco_recoil_E_calo,
    "reco_recoil_E_passive": reco_recoil_E_passive, "reco_anatool_E_nu": reco_anatool_E_nu, "reco_anatool_W": reco_anatool_W,
    "reco_visible_E": reco_visible_E, "reco_n_prongs": reco_n_prongs, "reco_n_hadron_tracks": reco_n_hadron_tracks,
    "reco_E_nu_calo": reco_E_nu_calo, "reco_q0_calo": reco_q0_calo, "reco_Q2_calo": reco_Q2_calo, "reco_q3_calo": reco_q3_calo,
    "reco_Q2_anatool": reco_Q2_anatool, "unit": unit,
    "reco_proton_p": reco_proton_p, "reco_proton_T": reco_proton_T, "reco_proton_theta": reco_proton_theta,
    "reco_proton_theta_deg": reco_proton_theta_deg, "reco_proton_pT": reco_proton_pT, "reco_proton_score1": reco_proton_score1,
    "reco_n_michel": reco_n_michel, "reco_n_iso_blobs": reco_n_iso_blobs,
    "reco_dpT": reco_dpT, "reco_dpTx": reco_dpTx, "reco_dpTy": reco_dpTy, "reco_dalphaT": reco_dalphaT,
    "reco_dalphaT_deg": reco_dalphaT_deg, "reco_dphiT": reco_dphiT, "reco_dphiT_deg": reco_dphiT_deg,
    "reco_dpL": reco_dpL, "reco_pn": reco_pn,
}


def namespace(r: dict, params: dict | None = None) -> dict:
    """Reco columns + reco observables, the latter evaluated lazily on first use by an expression."""
    from .observables import LazyNamespace
    ns = {k: np.asarray(v) for k, v in r.items() if k != "__meta__"}
    ns["M_MU"] = M_MU
    return LazyNamespace(ns, {name: (lambda f=f: f(r, params=params)) for name, f in RECO_OBSERVABLES.items()})


def evaluate(name_or_expr: str, r: dict, params: dict | None = None) -> np.ndarray:
    s = str(name_or_expr).strip()
    if s in RECO_OBSERVABLES:
        return RECO_OBSERVABLES[s](r, params=params)
    if is_identifier(s):
        return _col(r, s)
    n = len(next(iter(v for k, v in r.items() if k != "__meta__")))
    return eval_expr(s, namespace(r, params), n)
