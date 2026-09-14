"""Transverse kinematic imbalance observables (ndp.channels.observables: dpT, dpTx, dpTy, dalphaT,
dphiT, dpL, pn) on hand-built muon + proton events, and the p_n reconstruction identity."""
import numpy as np
from _helpers import site, have_mc_cache, skip
from ndp.channels import load_channel
from ndp.channels import observables as obs
from ndp.channels.observables import NUMI_BEAM_ANGLE_RAD
from ndp.events import TruthTable, M_P, M_MU, M_N

SPEC = {"type": "minerva_ccqelike_1mu1p", "nu_pdg": [14], "current": "CC", "lep_pdg": 13,
        "muon": {"theta_max_deg": 17.0, "p_min_gev": 2.0, "p_max_gev": 20.0},
        "proton": {"theta_max_deg": 70.0, "p_min_gev": 0.5, "p_max_gev": 1.1, "min_count": 1}}
PARAMS = {"tki": {"m_A_gev": 11.174864, "m_n_gev": 0.93956542, "excitation_b_gev": 0.02713}}
M_A = PARAMS["tki"]["m_A_gev"]
M_AP = M_A - PARAMS["tki"]["m_n_gev"] + PARAMS["tki"]["excitation_b_gev"]


def _table(events):
    """events = list of (muon 3-vector, [proton 3-vectors]) in GeV, detector frame, all numu CC."""
    n = len(events)
    mus = np.array([e[0] for e in events], float)
    parts = [(2212, p) for e in events for p in e[1]]
    fs_p = np.array([p for _, p in parts], float).reshape(-1, 3)
    cols = {
        "nu_pdg": np.full(n, 14), "E_nu": np.full(n, 6.0), "lep_pdg": np.full(n, 13),
        "lep_px": mus[:, 0], "lep_py": mus[:, 1], "lep_pz": mus[:, 2], "lep_E": np.sqrt((mus ** 2).sum(1) + M_MU ** 2),
        "current": np.ones(n, int), "int_type": np.ones(n, int), "target_Z": np.full(n, 6), "target_A": np.full(n, 12),
        "Q2": np.full(n, 0.3), "W": np.full(n, 1.0), "weight": np.ones(n),
        "fs_offsets": np.concatenate([[0], np.cumsum([len(e[1]) for e in events])]).astype(np.int64),
        "fs_pdg": np.array([pdg for pdg, _ in parts], int), "fs_E": np.sqrt((fs_p ** 2).sum(1) + M_P ** 2),
        "fs_px": fs_p[:, 0], "fs_py": fs_p[:, 1], "fs_pz": fs_p[:, 2],
    }
    return TruthTable(cols, {"source": "unit", "has_geometry": False})


def _ev(name, t, frame="detector"):
    return obs.evaluate(name, t, frame=frame, signal=SPEC, params=PARAMS)


def test_transverse_definitions_on_textbook_configurations():
    t = _table([
        ((0.0, -0.4, 4.0), [(0.0, 0.4, 0.6)]),        # 0 exactly back-to-back, equal pT
        ((0.0, -0.4, 4.0), [(0.3, 0.4, 0.6)]),        # 1 proton pT = -muon pT + 0.3 x-hat
        ((0.0, -0.5, 4.0), [(0.0, 0.3, 0.6)]),        # 2 proton decelerated (less pT, still back-to-back)
        ((0.0, -0.5, 4.0), [(0.0, 0.7, 0.6)]),        # 3 proton accelerated
        ((0.0, -0.4, 4.0), []),                       # 4 no proton
    ])
    dpt, dptx, dpty = _ev("dpT", t), _ev("dpTx", t), _ev("dpTy", t)
    da, dphi = _ev("dalphaT_deg", t), _ev("dphiT_deg", t)
    assert np.allclose([dpt[0], dptx[0], dpty[0], dphi[0]], 0.0, atol=1e-12) and np.isnan(da[0])
    assert np.allclose([dpt[1], dptx[1], dpty[1], da[1]], [0.3, 0.3, 0.0, 90.0]) and np.isclose(dphi[1], np.degrees(np.arccos(0.8)))
    assert np.allclose([dpt[2], dptx[2], dpty[2], da[2], dphi[2]], [0.2, 0.0, -0.2, 180.0, 0.0], atol=1e-9)
    assert np.allclose([dpt[3], dpty[3], da[3]], [0.2, 0.2, 0.0], atol=1e-9)
    assert all(np.isnan(v[4]) for v in (dpt, dptx, dpty, da, dphi, _ev("dpL", t), _ev("pn", t)))


def test_dptx_sign_is_z_cross_muon_hat():
    # muon along -y: z x (-y) = +x, so a proton pushed to +x gives dpTx > 0 (arXiv:2503.15047 Fig. 1);
    # muon along +x: z x (+x) = +y, so a proton pushed to +y gives dpTx > 0.
    t = _table([((0.0, -0.4, 4.0), [(0.3, 0.4, 0.6)]), ((0.4, 0.0, 4.0), [(-0.4, 0.3, 0.6)]), ((0.4, 0.0, 4.0), [(-0.4, -0.3, 0.6)])])
    assert np.allclose(_ev("dpTx", t), [0.3, 0.3, -0.3])


def test_pn_identity_for_one_nucleon_knockout_on_a_nucleus():
    """A QE event with initial neutron momentum p_n and remnant mass m_A' conserving 4-momentum:
    dpL must return p_n,z and pn must return |p_n| exactly (Furmanski & Sobczyk)."""
    p_n = np.array([0.10, 0.15, 0.05]); mu = np.array([0.1, -0.5, 3.5])   # solved: p_nu ~ 3.83 GeV, |p_p| ~ 0.75, theta_p ~ 60 deg
    E_mu = np.sqrt(mu @ mu + M_MU ** 2); E_rem = np.sqrt(M_AP ** 2 + p_n @ p_n)
    def resid(pnu):
        pp = np.array([0.0, 0.0, pnu]) + p_n - mu
        return pnu + M_A - E_mu - np.sqrt(pp @ pp + M_P ** 2) - E_rem
    lo, hi = 1.0, 20.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if resid(mid) < 0 else (lo, mid)
    pnu = 0.5 * (lo + hi)
    pp = np.array([0.0, 0.0, pnu]) + p_n - mu
    assert 0.5 < np.linalg.norm(pp) < 1.1 and np.degrees(np.arccos(pp[2] / np.linalg.norm(pp))) < 70
    t = _table([(tuple(mu), [tuple(pp)])])
    assert np.isclose(_ev("dpT", t)[0], np.hypot(p_n[0], p_n[1]))
    assert np.isclose(_ev("dpL", t)[0], p_n[2], atol=1e-9)
    assert np.isclose(_ev("pn", t)[0], np.linalg.norm(p_n), atol=1e-9)


def test_beam_frame_rotation_is_applied_consistently_to_muon_and_proton():
    a = NUMI_BEAM_ANGLE_RAD
    def unrotate(v):   # detector-frame vector whose beam-frame image is v
        x, y, z = v
        return (x, y * np.cos(a) + z * np.sin(a), -y * np.sin(a) + z * np.cos(a))
    mu, pp = (0.2, -0.4, 3.0), (-0.1, 0.55, 0.45)
    t_det = _table([(mu, [pp])])
    t_beam = _table([(unrotate(mu), [unrotate(pp)])])
    for name in ("dpT", "dpTx", "dpTy", "dalphaT", "dphiT", "dpL", "pn", "proton_p", "proton_theta"):
        assert np.isclose(_ev(name, t_det, "detector")[0], _ev(name, t_beam, "beam")[0], atol=1e-9), name
    # and the two frames genuinely differ on the same table
    assert not np.isclose(_ev("dpTy", t_det, "detector")[0], _ev("dpTy", t_det, "beam")[0])


def test_masses_come_from_the_manifest_not_the_code():
    t = _table([((0.0, -0.4, 4.0), [(0.3, 0.4, 0.6)])])
    assert np.isfinite(obs.evaluate("dpT", t, frame="detector", signal=SPEC)[0])
    try:
        obs.evaluate("dpL", t, frame="detector", signal=SPEC)
        raise AssertionError("dpL must refuse to run without observable_params.tki")
    except KeyError:
        pass
    ch = load_channel("minerva_me_ccqelike_1mu1p")
    assert ch.observable_params["tki"]["m_A_gev"] == M_A
    assert np.isfinite(ch.evaluate("pn", t)[0]) and np.isfinite(ch.evaluate("dalphaT_deg*1.0", t)[0])


def test_tki_on_the_open_data_mc_matches_the_signal_run():
    """Pins runs/2026-09-14_signal_minerva_me_ccqelike_1mu1p_2/summary.json (tki_fiducial_signal)."""
    if not have_mc_cache():
        skip("needs the cached MC tables")
    cfg = site(); ch = load_channel("minerva_me_ccqelike_1mu1p")
    t = TruthTable.load(cfg.data_dir / "cache/truth_mc110040.npz")
    sel = ch.is_signal(t) & ch.in_phase_space(t)
    assert int(sel.sum()) == 5392
    v = {k: ch.evaluate(k, t)[sel] for k in ("dpT", "dpTy", "dalphaT_deg", "pn")}
    assert all(np.isfinite(x).all() for x in v.values())
    assert (v["dalphaT_deg"] >= 0).all() and (v["dalphaT_deg"] <= 180).all()
    # medians measured by that run (summary.json: tki_fiducial_signal.<name>.median)
    ref = {"dpT": 0.26107235191911404, "dpTy": -0.1271141526271669, "dalphaT_deg": 132.90949287599983, "pn": 0.3457655167595715}
    for k, m in ref.items():
        assert np.isclose(np.median(v[k]), m, rtol=1e-9, atol=1e-9), (k, np.median(v[k]))
