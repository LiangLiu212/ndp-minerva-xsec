"""Reconstruction-level selections (ndp.channels.selections, adapters/minerva_anatuple.py::
ccqelike_1mu1p_cutflow): rule-by-rule on synthetic cache columns, registry behaviour, and the
certified counts of the 1mu1p selection on the open-data caches (runs/2026-09-14_selection_*_2)."""
import numpy as np
from _helpers import site, have_mc_cache, have_data_file, skip
from ndp.channels import load_channel
from ndp.channels.selections import select, cutflow, SELECTIONS
from ndp.adapters.minerva_anatuple import ccqelike_1mu1p_cutflow, CCQELIKE_1MU1P_LABELS, load_reco_cache, cache_tag
from ndp.events import TruthTable

PARAMS = {"muon": {"theta_max_deg": 17.0, "p_min_gev": 2.0, "p_max_gev": 20.0},
          "proton": {"theta_max_deg": 90.0, "p_min_gev": 0.4, "p_max_gev": 1.3},
          "proton_score1_min": 0.35, "require_contained": True, "n_michel_max": 0, "n_iso_blobs_max": 1}


def _table(n, **over):
    """A synthetic v3 cache in which every row passes; `over` sets a column to a list."""
    base = {"reco_vtx_x": 0.0, "reco_vtx_y": 0.0, "reco_vtx_z": 7000.0, "reco_minos_matched": True, "reco_n_dead_discr": 0.0,
            "reco_minos_qp": -0.2, "reco_theta": np.deg2rad(8.0), "reco_p": 5.0, "reco_proton_p": 0.8,
            "reco_proton_theta_beam": np.deg2rad(40.0), "reco_proton_score1": 0.9, "reco_proton_exiting": 0.0,
            "reco_n_michel": 0.0, "reco_n_iso_blobs": 0.0}
    r = {k: np.full(n, v) for k, v in base.items()}
    for k, v in over.items():
        r[k] = np.asarray(v)
    r["__meta__"] = {"cache_version": 3}
    return r


def test_each_cut_fails_exactly_where_intended():
    n = 13
    r = _table(n)
    r["reco_vtx_z"][1] = 5000.0                    # ZRange
    r["reco_vtx_x"][2] = 900.0                     # Apothem
    r["reco_minos_matched"][3] = False             # HasMINOSMatch
    r["reco_n_dead_discr"][4] = 2.0                # NoDeadtime
    r["reco_minos_qp"][5] = +0.2                   # IsNeutrino
    r["reco_theta"][6] = np.deg2rad(18.0)          # MuonWindow
    r["reco_proton_p"][7] = np.nan                 # HasProtonCandidate (no candidate -> NaN)
    r["reco_proton_exiting"][8] = 1.0              # ProtonContained
    r["reco_proton_score1"][9] = 0.2               # ProtonScore
    r["reco_proton_p"][10] = 1.4                   # ProtonWindow
    r["reco_n_michel"][11] = 1.0                   # NoMichel
    r["reco_n_iso_blobs"][12] = 2.0                # IsoBlobs
    passed, failing = ccqelike_1mu1p_cutflow(r, PARAMS)
    assert passed.tolist() == [True] + [False] * 12
    assert failing.tolist() == [-1] + list(range(12))
    assert len(CCQELIKE_1MU1P_LABELS) == 12
    # containment switched off in the manifest -> row 8 survives
    p2, _ = ccqelike_1mu1p_cutflow(r, dict(PARAMS, require_contained=False))
    assert p2[8] and p2.sum() == 2


def test_windows_and_thresholds_are_read_from_params():
    r = _table(4, reco_p=[1.9, 2.1, 19.9, 20.1], reco_proton_score1=[0.34, 0.36, 0.9, 0.9])
    p, f = ccqelike_1mu1p_cutflow(r, PARAMS)
    assert p.tolist() == [False, True, True, False] and f[0] == 5 and f[3] == 5
    r = _table(3, reco_proton_theta_beam=np.deg2rad([89.0, 91.0, 10.0]), reco_proton_p=[0.5, 0.5, 0.39])
    p, f = ccqelike_1mu1p_cutflow(r, PARAMS)
    assert p.tolist() == [True, False, False] and f[1] == 9 and f[2] == 9
    r = _table(2, reco_n_iso_blobs=[1.0, 2.0])
    assert ccqelike_1mu1p_cutflow(r, dict(PARAMS, n_iso_blobs_max=2))[0].all()


def test_registry_dispatch_and_inclusive_equivalence():
    ch = load_channel("minerva_me_cc_inclusive_ptpz")
    r = {"passed": np.array([True, False, True]), "failing_cut": np.array([-1, 2, -1]), "reco_p": np.ones(3)}
    assert np.array_equal(select(ch, r), r["passed"])
    steps = cutflow(ch, r)
    assert steps[0][0] == "all_candidates" and steps[0][1].all()
    assert [s for s, _ in steps][1:] == ["ZRange", "Apothem", "MaxMuonAngle", "HasMINOSMatch", "NoDeadtime", "IsNeutrino"]
    assert steps[3][1].tolist() == [True, False, True]      # failed at cut index 2 -> out from the third step on
    r2 = {"passed": np.array([True, False])}                 # toy tables carry only the flag
    assert np.array_equal(select(ch, r2), r2["passed"]) and cutflow(ch, r2)[-1][0] == "passed"
    ch1 = load_channel("minerva_me_ccqelike_1mu1p")
    assert ch1.selection["name"] in SELECTIONS
    try:
        select(ch1, {"passed": np.ones(2, bool), "__meta__": {"cache_version": 2}})
        raise AssertionError("a v2 cache must be refused")
    except ValueError:
        pass


def test_certified_counts_on_the_open_data_caches():
    """Pins runs/2026-09-14_selection_minerva_me_ccqelike_1mu1p_2 (cutflow.json, summary.json)."""
    if not (have_mc_cache() and have_data_file()):
        skip("needs the cached data + MC tables")
    cfg = site(); ch = load_channel("minerva_me_ccqelike_1mu1p"); cache = cfg.data_dir / "cache"
    rd = load_reco_cache(cache / f"reco_{cache_tag(ch.data['reco_data_files'][0])}.npz")
    rm = load_reco_cache(cache / f"reco_{cache_tag(ch.data['reco_mc_files'][0])}.npz")
    assert [int(m.sum()) for _, m in cutflow(ch, rd)] == [6304, 2658, 1656, 913, 899, 860, 766, 414, 402, 128, 123, 99, 64]
    assert [int(m.sum()) for _, m in cutflow(ch, rm)] == [186205, 126254, 80078, 48189, 47562, 44506, 40155, 22618, 21929, 8096, 7639, 6002, 3904]
    rt = TruthTable.load(cache / f"reco_{cache_tag(ch.data['reco_mc_files'][0])}_truthcols.npz")
    assert rt.has_fs
    sel = select(ch, rm); sig = ch.is_signal(rt) & ch.in_phase_space(rt)
    assert int((sel & sig).sum()) == 1897 and abs((sel & sig).sum() / sel.sum() - 0.4859118852459016) < 1e-12
    truth = TruthTable.load(cache / "truth_mc110040.npz")
    assert int((ch.is_signal(truth) & ch.in_phase_space(truth)).sum()) == 5392
    # the inclusive channel's certified selection is untouched by the v3 cache
    assert int(rm["passed"].sum()) == 43643 and int(rd["passed"].sum()) == 844
