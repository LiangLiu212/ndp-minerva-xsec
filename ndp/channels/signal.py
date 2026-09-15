"""Truth-level signal definitions: which generated events a channel counts as signal.

A channel manifest's `signal:` block names a definition with `type:` and carries its
parameters; the code here only evaluates them. Every parameter that is a physics choice
(kinematic windows, veto thresholds) therefore lives in the YAML, next to its `status`.

Registered types
----------------
cc_lepton (default when `type` is absent)
    {nu_pdg: [14], current: CC}  — the CC-inclusive definition the certified
    `minerva_me_cc_inclusive_ptpz` channel uses (unchanged behaviour).

minerva_ccqelike_1mu1p
    MINERvA's CCQE-like "muon + proton" definition (arXiv:2503.15047, Sec. "Analysis and
    results", verbatim in the MINERvA repo's `papers/minerva/2503.15047/paper_2503.15047.md`
    Sec. 3): one mu- inside a (theta, p) window, at least one proton inside a (theta, p)
    window, and no mesons, no baryons heavier than the neutron, no photons above a threshold.
    The highest-momentum proton inside the window is the "leading proton" that proton and TKI
    observables are built from. Parameters (channel YAML):

        signal:
          type: minerva_ccqelike_1mu1p
          nu_pdg: [14]
          current: CC
          lep_pdg: 13
          muon:   {theta_max_deg: 17.0, p_min_gev: 2.0, p_max_gev: 20.0}
          proton: {theta_max_deg: 70.0, p_min_gev: 0.5, p_max_gev: 1.1, min_count: 1}
          veto:   {mesons: true, heavy_baryons: true, photon_E_max_gev: 0.010}

    Angles are measured in the channel's `phase_space.frame` (the paper: "with respect to the
    beam"), applied to the final-state protons with the same rotation as the lepton.

Final-state particle classes (`fs_classes`)
------------------------------------------
    meson          100 <= |pdg| < 1000  (pi, K, eta, rho, K*, D, ...; the paper: "no mesons")
    heavy_baryon   1000 <= |pdg| < 10000 and |pdg| not in {2212, 2112}
                   (hyperons, charmed baryons, any undecayed resonance; antinucleons are *not*
                   heavier than the neutron and are not vetoed — same as MAT's `IsQELike`)
    photon         pdg == 22
    proton/neutron pdg == 2212 / 2112
    nucleus        |pdg| >= 1e9 (nuclear remnants) — ignored
    pseudo         |pdg| >= 2e9 (GENIE bookkeeping particles, e.g. 2000000101) — ignored
    charged_lepton, neutrino — ignored by the vetoes (the primary lepton is `lep_*`)

The MAT reference (`MAT-MINERvA/universes/CCQE3DFitsSystematics.cxx::IsQELike`) enumerates
PDG codes explicitly; on the open-data MC every code it lists falls in the classes above, so
the two agree on that sample, and the class rules also cover codes MAT's list omits.
"""
from __future__ import annotations

import numpy as np

from ..events import TruthTable, M_P
from .observables import NUMI_BEAM_ANGLE_RAD, frame_rotation_angle, lep_p, lep_theta, native_frame, rotate_about_x

NUCLEUS_MIN = 1_000_000_000
PSEUDO_MIN = 2_000_000_000


# ---- kinematics helpers ---------------------------------------------------------------------
def rotate_to_frame(px, py, pz, frame: str, native: str = "detector"):
    """Momentum components stored in the `native` frame -> the requested frame (same convention as
    the lepton: rotation about x by NUMI_BEAM_ANGLE_RAD for detector -> beam, the inverse back)."""
    return rotate_about_x(px, py, pz, frame_rotation_angle(native, frame))


def _theta(px, py, pz):
    p = np.sqrt(px * px + py * py + pz * pz)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.where(p > 0, pz / p, 1.0)
    return np.arccos(np.clip(c, -1.0, 1.0)), p


# ---- final-state classes ------------------------------------------------------------------
def fs_classes(t: TruthTable) -> dict:
    """Per-final-state-particle boolean arrays (see module docstring)."""
    if not t.has_fs:
        raise ValueError("signal definition needs final-state particles (fs_* columns)")
    pdg = t["fs_pdg"]
    apdg = np.abs(pdg)
    nucleus = (apdg >= NUCLEUS_MIN) & (apdg < PSEUDO_MIN)
    pseudo = apdg >= PSEUDO_MIN
    baryon = (apdg >= 1000) & (apdg < 10000)
    return {
        "meson": (apdg >= 100) & (apdg < 1000),
        "heavy_baryon": baryon & ~np.isin(apdg, [2212, 2112]),
        "photon": pdg == 22,
        "proton": pdg == 2212,
        "neutron": pdg == 2112,
        "charged_lepton": np.isin(apdg, [11, 13, 15]),
        "neutrino": np.isin(apdg, [12, 14, 16]),
        "nucleus": nucleus,
        "pseudo": pseudo,
    }


def fs_count(t: TruthTable, particle_mask: np.ndarray) -> np.ndarray:
    """Per-event number of final-state particles satisfying `particle_mask`."""
    return t.fs_sum(np.ones(len(particle_mask)), particle_mask).astype(np.int64)


# ---- derived columns (tables whose final-state list was replaced by per-event summaries) ----
#: per-event columns a skim stores instead of fs_* (see adapters.minerva_anatuple.derive_fs_columns)
DERIVED_LP = ("lp_p", "lp_theta", "lp_pT", "lp_px", "lp_py", "lp_pz", "lp_E", "lp_n_in_window")
DERIVED_VETO = ("n_meson", "n_heavy_baryon", "n_photon_hard")
_WINDOW_KEYS = ("theta_max_deg", "p_min_gev", "p_max_gev")


def _same_window(a: dict | None, b: dict | None) -> bool:
    a, b = a or {}, b or {}
    return all((a.get(k) is None) == (b.get(k) is None) and (a.get(k) is None or float(a[k]) == float(b[k])) for k in _WINDOW_KEYS)


def derived_info(t: TruthTable) -> dict | None:
    """The `derived` block of a skim's meta (frame, proton_window, photon_E_max_gev, channel), or None."""
    d = t.meta.get("derived")
    return dict(d) if isinstance(d, dict) else None


def _require_derived(t: TruthTable, columns, frame: str | None = None, window: dict | None = None,
                     photon_E_max_gev: float | None = None) -> dict:
    d = derived_info(t)
    missing = [c for c in columns if c not in t]
    if d is None or missing:
        raise ValueError("table has no fs_* columns and no usable derived columns "
                         f"(missing {missing or 'the meta.derived block'}); build the skim with derive_fs_columns")
    if frame is not None and d.get("frame") != frame:
        raise ValueError(f"derived columns were computed in frame {d.get('frame')!r}, requested {frame!r}")
    if window is not None and not _same_window(d.get("proton_window"), window):
        raise ValueError(f"derived columns use proton window {d.get('proton_window')}, requested {window}")
    if photon_E_max_gev is not None and float(d.get("photon_E_max_gev", -1)) != float(photon_E_max_gev):
        raise ValueError(f"derived photon threshold {d.get('photon_E_max_gev')} != requested {photon_E_max_gev}")
    return d


def veto_counts(t: TruthTable, photon_E_max_gev: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-event (n_meson, n_heavy_baryon, n_photon_above_threshold) from fs_* or from derived columns."""
    if t.has_fs:
        cls = fs_classes(t)
        hard = cls["photon"] & (t["fs_E"] > float(photon_E_max_gev))
        return fs_count(t, cls["meson"]), fs_count(t, cls["heavy_baryon"]), fs_count(t, hard)
    _require_derived(t, DERIVED_VETO, photon_E_max_gev=photon_E_max_gev)
    return (t["n_meson"].astype(np.int64), t["n_heavy_baryon"].astype(np.int64), t["n_photon_hard"].astype(np.int64))


# ---- leading proton -------------------------------------------------------------------------
def leading_proton(t: TruthTable, frame: str = "detector", window: dict | None = None) -> dict:
    """Highest-momentum final-state proton inside `window` per event.

    window = {theta_max_deg, p_min_gev, p_max_gev} (any key may be absent). Returns per-event
    arrays: `index` (flat fs index, -1 if none), `n_in_window`, `p`, `theta` [rad], `pT`,
    `px`, `py`, `pz`, `E` (GeV, in `frame`; NaN where the event has no proton in the window).
    On a skim (no fs_* columns) the stored `lp_*` columns are returned, provided they were
    derived in the same frame and window; `index` is then 0 / -1 (has / has no proton).
    """
    if not t.has_fs:
        _require_derived(t, DERIVED_LP, frame=frame, window=window)
        n_in = t["lp_n_in_window"].astype(np.int64)
        return {"index": np.where(n_in > 0, 0, -1).astype(np.int64), "n_in_window": n_in, "p": t["lp_p"], "theta": t["lp_theta"],
                "pT": t["lp_pT"], "px": t["lp_px"], "py": t["lp_py"], "pz": t["lp_pz"], "E": t["lp_E"]}
    w = window or {}
    px, py, pz = rotate_to_frame(t["fs_px"], t["fs_py"], t["fs_pz"], frame, native_frame(t))
    theta, p = _theta(px, py, pz)
    ok = t["fs_pdg"] == 2212
    if "theta_max_deg" in w:
        ok &= theta < np.deg2rad(float(w["theta_max_deg"]))
    if "p_min_gev" in w:
        ok &= p > float(w["p_min_gev"])
    if "p_max_gev" in w:
        ok &= p < float(w["p_max_gev"])
    ev = t.fs_event_index()
    n_in = np.zeros(t.n, dtype=np.int64)
    np.add.at(n_in, ev[ok], 1)
    index = np.full(t.n, -1, dtype=np.int64)
    if ok.any():
        cand = np.flatnonzero(ok)
        order = np.lexsort((-p[cand], ev[cand]))          # by event, then descending momentum
        cand = cand[order]
        first = np.unique(ev[cand], return_index=True)[1]  # first row per event = leading proton
        index[ev[cand[first]]] = cand[first]
    has = index >= 0
    safe = np.where(has, index, 0)
    def pick(a):
        return np.where(has, a[safe], np.nan)
    return {"index": index, "n_in_window": n_in, "p": pick(p), "theta": pick(theta),
            "pT": pick(np.sqrt(px * px + py * py)), "px": pick(px), "py": pick(py), "pz": pick(pz),
            "E": pick(t["fs_E"]) if "fs_E" in t else np.where(has, np.sqrt(p[safe] ** 2 + M_P ** 2), np.nan)}


# ---- definitions ---------------------------------------------------------------------------
def _current_code(spec: dict) -> int:
    return {"CC": 1, "NC": 2}.get(str(spec.get("current", "CC")).upper(), 1)


def cc_lepton(spec: dict, t: TruthTable, frame: str = "detector") -> np.ndarray:
    nu = np.isin(t["nu_pdg"], spec.get("nu_pdg", [14]))
    return nu & (t["current"] == _current_code(spec))


def cutflow_1mu1p(spec: dict, t: TruthTable, frame: str = "detector") -> list[tuple[str, np.ndarray]]:
    """Ordered (label, cumulative mask) pairs of the minerva_ccqelike_1mu1p definition."""
    mu = spec.get("muon", {})
    pr = spec.get("proton", {})
    veto = spec.get("veto", {})
    steps: list[tuple[str, np.ndarray]] = []
    m = cc_lepton(spec, t, frame) & (t["lep_pdg"] == int(spec.get("lep_pdg", 13)))
    steps.append(("numu_cc_muon", m.copy()))
    th, p = lep_theta(t, frame), lep_p(t, frame)
    if "theta_max_deg" in mu:
        m &= th < np.deg2rad(float(mu["theta_max_deg"]))
    if "p_min_gev" in mu:
        m &= p > float(mu["p_min_gev"])
    if "p_max_gev" in mu:
        m &= p < float(mu["p_max_gev"])
    steps.append(("muon_window", m.copy()))
    lp = leading_proton(t, frame, pr)
    m &= lp["n_in_window"] >= int(pr.get("min_count", 1))
    steps.append(("proton_in_window", m.copy()))
    n_meson, n_heavy, n_hard = veto_counts(t, float(veto.get("photon_E_max_gev", np.inf)))
    if veto.get("mesons", True):
        m &= n_meson == 0
        steps.append(("no_mesons", m.copy()))
    if veto.get("heavy_baryons", True):
        m &= n_heavy == 0
        steps.append(("no_heavy_baryons", m.copy()))
    if "photon_E_max_gev" in veto:
        m &= n_hard == 0
        steps.append(("no_photons_above_threshold", m.copy()))
    return steps


def minerva_ccqelike_1mu1p(spec: dict, t: TruthTable, frame: str = "detector") -> np.ndarray:
    return cutflow_1mu1p(spec, t, frame)[-1][1]


SIGNAL_TYPES = {"cc_lepton": cc_lepton, "minerva_ccqelike_1mu1p": minerva_ccqelike_1mu1p}
CUTFLOWS = {"cc_lepton": lambda s, t, f: [("cc_lepton", cc_lepton(s, t, f))], "minerva_ccqelike_1mu1p": cutflow_1mu1p}


def signal_type(spec: dict) -> str:
    kind = str(spec.get("type", "cc_lepton"))
    if kind not in SIGNAL_TYPES:
        raise KeyError(f"unknown signal type {kind!r}; known: {sorted(SIGNAL_TYPES)}")
    return kind


def evaluate_signal(spec: dict, t: TruthTable, frame: str = "detector") -> np.ndarray:
    return SIGNAL_TYPES[signal_type(spec)](spec, t, frame)


def signal_cutflow(spec: dict, t: TruthTable, frame: str = "detector") -> list[tuple[str, np.ndarray]]:
    return CUTFLOWS[signal_type(spec)](spec, t, frame)
