"""Truth-level signal definitions (ndp.channels.signal): rule-by-rule on hand-built events, and
certification of the minerva_ccqelike_1mu1p counts on the cached open-data MC."""
import numpy as np
from _helpers import toy_truth, site, have_mc_cache, skip
from ndp.channels import load_channel
from ndp.channels import signal as sig
from ndp.channels import observables as obs
from ndp.events import TruthTable, M_P

SPEC = {"type": "minerva_ccqelike_1mu1p", "nu_pdg": [14], "current": "CC", "lep_pdg": 13,
        "muon": {"theta_max_deg": 17.0, "p_min_gev": 2.0, "p_max_gev": 20.0},
        "proton": {"theta_max_deg": 70.0, "p_min_gev": 0.5, "p_max_gev": 1.1, "min_count": 1},
        "veto": {"mesons": True, "heavy_baryons": True, "photon_E_max_gev": 0.010}}


def _p4(pdg, p, theta_deg, phi_deg=0.0):
    m = {2212: M_P, 2112: 0.93957, 211: 0.13957, -211: 0.13957, 111: 0.13498, 22: 0.0, 3122: 1.11568,
         -2212: M_P, 13: 0.10566, -13: 0.10566, 2000000101: 0.0, 1000060120: 11.17}.get(pdg, 0.5)
    th, ph = np.deg2rad(theta_deg), np.deg2rad(phi_deg)
    return pdg, np.sqrt(p * p + m * m), p * np.sin(th) * np.cos(ph), p * np.sin(th) * np.sin(ph), p * np.cos(th)


def _events(fs_lists, muon=(5.0, 5.0)):
    """Detector-frame events, each with a mu- of (p [GeV], theta [deg]) and the given FS particles."""
    n = len(fs_lists)
    p, th = muon
    parts = [q for ev in fs_lists for q in ev]
    cols = {
        "nu_pdg": np.full(n, 14), "E_nu": np.full(n, 6.0), "lep_pdg": np.full(n, 13),
        "lep_px": np.full(n, p * np.sin(np.deg2rad(th))), "lep_py": np.zeros(n), "lep_pz": np.full(n, p * np.cos(np.deg2rad(th))),
        "lep_E": np.full(n, np.sqrt(p * p + 0.10566 ** 2)), "current": np.ones(n, int), "int_type": np.ones(n, int),
        "target_Z": np.full(n, 6), "target_A": np.full(n, 12), "Q2": np.full(n, 0.3), "W": np.full(n, 1.0), "weight": np.ones(n),
        "fs_offsets": np.concatenate([[0], np.cumsum([len(ev) for ev in fs_lists])]).astype(np.int64),
        "fs_pdg": np.array([q[0] for q in parts], int), "fs_E": np.array([q[1] for q in parts]),
        "fs_px": np.array([q[2] for q in parts]), "fs_py": np.array([q[3] for q in parts]), "fs_pz": np.array([q[4] for q in parts]),
    }
    return TruthTable(cols, {"source": "unit", "has_geometry": False})


def test_rules_one_by_one_detector_frame():
    good = [_p4(13, 5.0, 5.0), _p4(2212, 0.8, 30.0), _p4(2112, 0.3, 80.0), _p4(22, 0.005, 10.0), _p4(2000000101, 0.01, 0.0)]
    fs = [
        good,                                             # 0 signal
        good + [_p4(211, 0.2, 20.0)],                     # 1 charged pion -> vetoed
        good + [_p4(111, 0.2, 20.0)],                     # 2 pi0 -> vetoed
        good + [_p4(3122, 0.9, 20.0)],                    # 3 Lambda -> vetoed
        good + [_p4(22, 0.020, 20.0)],                    # 4 20 MeV photon -> vetoed
        good + [_p4(-2212, 0.6, 20.0)],                   # 5 antiproton: not heavier than the neutron -> signal
        [_p4(13, 5.0, 5.0), _p4(2212, 0.4, 30.0)],        # 6 proton below 0.5 GeV/c -> no proton in window
        [_p4(13, 5.0, 5.0), _p4(2212, 1.2, 30.0)],        # 7 proton above 1.1 GeV/c
        [_p4(13, 5.0, 5.0), _p4(2212, 0.8, 75.0)],        # 8 proton beyond 70 deg
        [_p4(13, 5.0, 5.0), _p4(2212, 0.8, 30.0), _p4(2212, 1.5, 10.0)],   # 9 one in window + one above: signal (out-of-window protons do not veto)
        good + [_p4(1000060120, 0.05, 0.0)],              # 10 nuclear remnant ignored
    ]
    t = _events(fs)
    m = sig.evaluate_signal(SPEC, t, "detector")
    assert m.tolist() == [True, False, False, False, False, True, False, False, False, True, True]
    lp = sig.leading_proton(t, "detector", SPEC["proton"])
    assert np.isclose(lp["p"][9], 0.8) and lp["n_in_window"][9] == 1 and lp["index"][6] == -1 and np.isnan(lp["p"][6])
    steps = dict(sig.cutflow_1mu1p(SPEC, t, "detector"))
    assert steps["proton_in_window"].sum() == 8 and steps["no_mesons"].sum() == 6 and steps["no_heavy_baryons"].sum() == 5


def test_muon_window_and_lepton_flavour():
    fs = [[_p4(13, 5.0, 5.0), _p4(2212, 0.8, 30.0)]] * 4
    for muon, expect in [((5.0, 5.0), True), ((5.0, 18.0), False), ((1.5, 5.0), False), ((25.0, 5.0), False)]:
        t = _events(fs, muon=muon)
        assert bool(sig.evaluate_signal(SPEC, t, "detector")[0]) is expect, muon
    t = _events(fs); t.columns["nu_pdg"][:] = -14; t.columns["lep_pdg"][:] = -13
    assert not sig.evaluate_signal(SPEC, t, "detector").any()
    t = _events(fs); t.columns["current"][:] = 2
    assert not sig.evaluate_signal(SPEC, t, "detector").any()


def test_leading_proton_is_highest_momentum_in_window_and_frame_rotates_protons_like_the_lepton():
    fs = [[_p4(13, 5.0, 5.0), _p4(2212, 0.6, 20.0, 0.0), _p4(2212, 1.0, 60.0, 90.0), _p4(2212, 1.05, 69.9, 180.0)]]
    t = _events(fs)
    lp = sig.leading_proton(t, "detector", SPEC["proton"])
    assert np.isclose(lp["p"][0], 1.05) and lp["n_in_window"][0] == 3
    # in the beam frame the third proton (phi = 180, i.e. along -x... theta only) keeps |p| but its angle changes
    lpb = sig.leading_proton(t, "beam", SPEC["proton"])
    assert np.isclose(lpb["p"][0], 1.05) or np.isclose(lpb["p"][0], 1.0)   # |p| invariant, window membership may change
    # the rotation applied to protons is the lepton's rotation
    px, py, pz = sig.rotate_to_frame(t["lep_px"], t["lep_py"], t["lep_pz"], "beam")
    assert np.allclose(np.sqrt(px * px + py * py), obs.lep_pT(t, "beam")) and np.allclose(pz, obs.lep_pz(t, "beam"))


def test_default_type_is_the_cc_inclusive_definition():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    t = toy_truth(3000)
    assert sig.signal_type(ch.signal) == "cc_lepton"
    expect = np.isin(t["nu_pdg"], [14]) & (t["current"] == 1)
    assert np.array_equal(ch.is_signal(t), expect)
    assert [n for n, _ in ch.signal_cutflow(t)] == ["cc_lepton"]


def test_proton_observables_follow_the_channel_signal_block():
    ch = load_channel("minerva_me_ccqelike_1mu1p")
    assert sig.signal_type(ch.signal) == "minerva_ccqelike_1mu1p" and ch.frame == "beam"
    t = _events([[_p4(13, 5.0, 5.0), _p4(2212, 0.8, 30.0), _p4(2212, 1.5, 10.0)]])
    x, _ = ch.observables(t)                       # channel binning x = proton_p (leading proton inside the window)
    assert np.isclose(x[0], 0.8)
    no_window = obs.evaluate("proton_p", t, frame="beam")            # without a signal block: leading proton overall
    assert np.isclose(no_window[0], 1.5)
    expr = obs.evaluate("proton_p*cos(proton_theta)", t, frame="beam", signal=ch.signal)
    assert np.isfinite(expr[0]) and expr[0] < 0.8


def test_certified_counts_on_the_open_data_mc():
    """Pins runs/2026-09-14_signal_minerva_me_ccqelike_1mu1p (cutflow.json / summary.json)."""
    if not have_mc_cache():
        skip("needs the cached MC tables")
    cfg = site(); ch = load_channel("minerva_me_ccqelike_1mu1p")
    t = TruthTable.load(cfg.data_dir / "cache/truth_mc110040.npz")
    steps = dict(ch.signal_cutflow(t))
    assert [int(steps[k].sum()) for k in ("numu_cc_muon", "muon_window", "proton_in_window", "no_mesons",
                                           "no_heavy_baryons", "no_photons_above_threshold")] == \
        [397604, 224629, 74184, 21974, 21973, 21369]
    fid = ch.in_phase_space(t)
    assert int((steps["no_photons_above_threshold"] & fid).sum()) == 5392
    from ndp.diagnostics import mat_isqelike
    cls = sig.fs_classes(t)
    hard = cls["photon"] & (t["fs_E"] > 0.010)
    ours = (sig.fs_count(t, cls["meson"]) == 0) & (sig.fs_count(t, cls["heavy_baryon"]) == 0) & (sig.fs_count(t, hard) == 0)
    numu = steps["numu_cc_muon"]
    assert int(((ours != mat_isqelike(t)) & numu).sum()) == 6      # anti-Lambda / anti-K0 / extra FS muon events, see docs/decisions.md
