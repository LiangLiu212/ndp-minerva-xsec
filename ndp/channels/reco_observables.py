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


def reco_p(r): return _col(r, "reco_p")
def reco_theta(r): return _col(r, "reco_theta")
def reco_theta_deg(r): return np.rad2deg(_col(r, "reco_theta"))
def reco_pT(r): return _col(r, "reco_pT")
def reco_pz(r): return _col(r, "reco_pz")
def reco_E_mu(r): return _col(r, "reco_E_mu")
def reco_recoil_E(r): return _col(r, "reco_recoil_E")
def reco_recoil_E_calo(r): return _col(r, "reco_recoil_E_calo")
def reco_recoil_E_passive(r): return _col(r, "reco_recoil_E_passive")
def reco_anatool_E_nu(r): return _col(r, "reco_anatool_E_nu")
def reco_anatool_W(r): return _col(r, "reco_anatool_W")
def reco_visible_E(r): return _col(r, "reco_visible_E")
def reco_n_prongs(r): return _col(r, "reco_n_prongs")
def reco_n_hadron_tracks(r): return _col(r, "reco_n_hadron_tracks")


def reco_E_nu_calo(r):
    """Calorimetric neutrino energy: muon energy + MasterAnaDev_recoil_E (the CC-inclusive recoil)."""
    return reco_E_mu(r) + reco_recoil_E(r)


def reco_q0_calo(r):
    return reco_recoil_E(r)


def reco_Q2_calo(r):
    """Q^2 = 2 E_nu (E_mu - p_mu cos theta) - m_mu^2 with the calorimetric E_nu [GeV^2]."""
    return 2.0 * reco_E_nu_calo(r) * (reco_E_mu(r) - reco_p(r) * np.cos(reco_theta(r))) - M_MU ** 2


def reco_q3_calo(r):
    q0 = reco_q0_calo(r)
    return np.sqrt(np.maximum(reco_Q2_calo(r) + q0 * q0, 0.0))


def reco_Q2_anatool(r):
    """Q^2 from the analysis tool's own E_nu (MasterAnaDev_E) and the muon [GeV^2]."""
    return 2.0 * reco_anatool_E_nu(r) * (reco_E_mu(r) - reco_p(r) * np.cos(reco_theta(r))) - M_MU ** 2


def unit(r):
    n = len(next(iter(r.values())))
    return np.full(n, 0.5)


RECO_OBSERVABLES = {
    "reco_p": reco_p, "reco_theta": reco_theta, "reco_theta_deg": reco_theta_deg, "reco_pT": reco_pT, "reco_pz": reco_pz,
    "reco_E_mu": reco_E_mu, "reco_recoil_E": reco_recoil_E, "reco_recoil_E_calo": reco_recoil_E_calo,
    "reco_recoil_E_passive": reco_recoil_E_passive, "reco_anatool_E_nu": reco_anatool_E_nu, "reco_anatool_W": reco_anatool_W,
    "reco_visible_E": reco_visible_E, "reco_n_prongs": reco_n_prongs, "reco_n_hadron_tracks": reco_n_hadron_tracks,
    "reco_E_nu_calo": reco_E_nu_calo, "reco_q0_calo": reco_q0_calo, "reco_Q2_calo": reco_Q2_calo, "reco_q3_calo": reco_q3_calo,
    "reco_Q2_anatool": reco_Q2_anatool, "unit": unit,
}


def namespace(r: dict) -> dict:
    ns = {k: np.asarray(v) for k, v in r.items() if k != "__meta__"}
    for name, f in RECO_OBSERVABLES.items():
        try:
            ns[name] = f(r)
        except KeyError:
            pass
    ns["M_MU"] = M_MU
    return ns


def evaluate(name_or_expr: str, r: dict) -> np.ndarray:
    s = str(name_or_expr).strip()
    if s in RECO_OBSERVABLES:
        return RECO_OBSERVABLES[s](r)
    if is_identifier(s):
        return _col(r, s)
    n = len(next(iter(v for k, v in r.items() if k != "__meta__")))
    return eval_expr(s, namespace(r), n)
