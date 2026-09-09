"""Measurements: user-defined observable pairs, truth/reco expressions, 1D grids, surrogate closure."""
import tempfile
from pathlib import Path

import numpy as np
from _helpers import toy_truth, site, have_mc_cache, skip
from ndp.channels import load_channel, load_measurement, list_measurements, published_measurement
from ndp.channels import observables as obs
from ndp.channels import reco_observables as robs
from ndp.channels.binning import Binning
from ndp.surrogate.binned import BinnedResponse
from ndp.surrogate.build import training_arrays, closure
from ndp.events import TruthTable


def test_published_measurement_is_the_channel_grid():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    m = published_measurement(ch)
    assert m.binning == ch.binning and m.release == "2106.16210" and m.name == "published"
    assert (m.x.truth, m.x.reco, m.y.truth, m.y.reco) == ("lep_pT", "reco_pT", "lep_pz", "reco_pz")
    assert "published" in list_measurements(ch)


def test_example_measurements_load():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    names = list_measurements(ch)
    for n in ("muon_p_theta", "enu_calorimetric", "q2_calorimetric"):
        assert n in names
        m = load_measurement(ch, n)
        assert m.release is None and m.status == "example"
    m1 = load_measurement(ch, "enu_calorimetric")
    assert m1.is_1d and m1.binning.n_y == 1 and m1.y.truth == "unit"
    m2 = load_measurement(ch, "muon_p_theta")
    assert not m2.is_1d and m2.binning.n_cells == 14 * 9


def test_truth_expressions_and_registry():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    t = toy_truth(3000)
    frame = ch.phase_space["frame"]
    assert np.allclose(obs.evaluate("sqrt(Q2)", t, frame=frame), np.sqrt(t["Q2"]))
    assert np.allclose(obs.evaluate("E_nu - lep_E", t, frame=frame), obs.q0(t))
    assert np.allclose(obs.evaluate("lep_p*cos(lep_theta)", t, frame=frame), obs.lep_pz(t, frame))
    assert np.allclose(obs.evaluate("unit", t), 0.5)
    try:
        obs.evaluate("no_such_thing + 1", t)
        assert False, "unknown name must raise"
    except KeyError:
        pass


def test_reco_expressions_and_registry():
    n = 500
    rng = np.random.default_rng(1)
    p = rng.uniform(1.5, 20, n); th = rng.uniform(0, 0.35, n)
    r = {"reco_p": p, "reco_theta": th, "reco_pT": p * np.sin(th), "reco_pz": p * np.cos(th),
         "reco_E_mu": np.sqrt(p ** 2 + 0.1057 ** 2), "reco_recoil_E": rng.exponential(1.0, n)}
    assert np.allclose(robs.evaluate("reco_E_nu_calo", r), r["reco_E_mu"] + r["reco_recoil_E"])
    assert np.allclose(robs.evaluate("reco_E_mu + reco_recoil_E", r), r["reco_E_mu"] + r["reco_recoil_E"])
    assert np.allclose(robs.evaluate("rad2deg(reco_theta)", r), np.rad2deg(th))
    q2 = robs.evaluate("reco_Q2_calo", r)
    assert np.all(q2 > -1e-9) and np.allclose(q2, 2 * (r["reco_E_mu"] + r["reco_recoil_E"]) * (r["reco_E_mu"] - p * np.cos(th)) - 0.1056583755 ** 2)
    try:
        robs.evaluate("reco_visible_E", r)
        assert False, "missing cache column must raise"
    except KeyError as e:
        assert "ndp data cache" in str(e)


def test_measurement_from_yaml_path_and_1d_binning():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    d = tempfile.mkdtemp()
    p = Path(d) / "my_obs.yaml"
    p.write_text("name: my_obs\nchannel: minerva_me_cc_inclusive_ptpz\n"
                 "x: {observable: 'E_nu - lep_E', reco: 'reco_recoil_E', units: GeV, edges: [0, 0.5, 1, 2, 5, 20]}\n")
    m = load_measurement(ch, p)
    assert m.name == "my_obs" and m.is_1d and m.binning.n_cells == 5 and m.status == "user"
    t = toy_truth(4000)
    x, y = m.truth_observables(ch, t)
    assert np.allclose(x, obs.q0(t)) and np.all(y == 0.5)
    sumw, sumw2, n_out, mask = m.truth_cells(ch, t)
    assert sumw.sum() + n_out == mask.sum()
    # a measurement belonging to another channel is refused
    p2 = Path(d) / "wrong.yaml"
    p2.write_text("name: wrong\nchannel: other\nx: {observable: E_nu, reco: reco_p, edges: [0, 1]}\n")
    try:
        load_measurement(ch, p2); assert False
    except ValueError:
        pass


def test_surrogate_closure_on_a_user_measurement_toy():
    """training_arrays + BinnedResponse.fit + closure on a toy paired sample, through the Measurement API."""
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    d = tempfile.mkdtemp(); p = Path(d) / "pth.yaml"
    p.write_text("name: pth\nchannel: minerva_me_cc_inclusive_ptpz\n"
                 "x: {observable: lep_p, reco: reco_p, edges: [1.5, 3, 5, 8, 12, 20, 40]}\n"
                 "y: {observable: lep_theta_deg, reco: 'rad2deg(reco_theta)', edges: [0, 5, 10, 15, 20]}\n")
    m = load_measurement(ch, p)
    truth = toy_truth(30000, seed=7)
    truth.meta["has_geometry"] = True
    rng = np.random.default_rng(3)
    # 'reco rows': a subset of the truth with smeared reco quantities (all rows selected, all signal)
    idx = np.sort(rng.choice(truth.n, 18000, replace=False))
    rt = truth.select(np.isin(np.arange(truth.n), idx))
    frame = ch.phase_space["frame"]
    p_true, th_true = obs.lep_p(rt, frame), obs.lep_theta(rt, frame)
    reco = {"passed": np.ones(rt.n, bool), "reco_p": p_true * rng.normal(1, 0.08, rt.n), "reco_theta": np.abs(th_true + rng.normal(0, 0.01, rt.n))}
    arrs = training_arrays(ch, m, truth, reco, rt)
    sur = BinnedResponse.fit(m.binning, **arrs["kw"], **arrs["bkg"], pot_mc=1.0)
    truth.meta["norm"] = {"kind": "pot", "pot": 1.0}
    clo = closure(ch, m, sur, truth, reco)
    assert clo["exact"], clo
    assert abs(clo["n_pred"] - clo["n_reco"]) < 1e-6


def test_user_measurement_surrogates_close_on_the_minerva_mc():
    """Every built surrogate under surrogates/<channel>/<measurement>/binned_* reproduces the training reco counts."""
    if not have_mc_cache():
        skip("needs the cached MC tables")
    from ndp.surrogate.base import load_surrogate
    from ndp.surrogate.build import load_training_cache
    cfg = site(); ch = load_channel("minerva_me_cc_inclusive_ptpz")
    try:
        tc = load_training_cache(cfg, ch)
    except ValueError as e:
        skip(str(e))
    checked = 0
    for name in list_measurements(ch):
        m = load_measurement(ch, name)
        for sdir in sorted(m.surrogate_root(cfg).glob("binned_mc*")):
            s = load_surrogate(sdir)
            clo = closure(ch, m, s, tc["truth"], tc["reco"])
            assert clo["exact"], (name, sdir, clo)
            checked += 1
    if checked == 0:
        skip("no measurement surrogates built")
