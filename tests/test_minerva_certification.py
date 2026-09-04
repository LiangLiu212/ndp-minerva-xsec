"""Certification against the MINERvA exploration repo's audited numbers (skipped without data)."""
import numpy as np
from _helpers import site, have_minerva_repo, have_mc_cache, have_data_file, skip
from ndp.channels import load_channel
from ndp.events import TruthTable


def test_selection_parity_vs_tool():
    if not (have_minerva_repo() and have_data_file()):
        skip("needs the MINERvA repo and the data AnaTuple")
    from ndp.adapters.minerva_anatuple import parity_vs_tool
    cfg = site()
    r = parity_vs_tool(cfg.data_dir / "MasterAnaDev_data_AnaTuple_run00010066_Playlist.root", cfg.minerva_repo, n=6304)
    assert r["n_mismatch"] == 0 and r["n_passed"] == 844


def test_audited_run_counts_detector_frame():
    """runs/2026-06-19_me_inclusive_ddxsec used the detector frame for truth; reproduce its counts."""
    if not have_mc_cache():
        skip("needs the cached MC tables")
    cfg = site(); ch = load_channel("minerva_me_cc_inclusive_ptpz")
    ch.phase_space = dict(ch.phase_space, frame="detector")
    truth = TruthTable.load(cfg.data_dir / "cache/truth_mc110040.npz")
    sig = ch.is_signal(truth) & ch.in_phase_space(truth)
    x, y = ch.observables(truth); g = ch.binning.digitize(x, y)
    assert int(sig.sum()) == 65041 and int((sig & (g >= 0)).sum()) == 65003
    z = np.load(cfg.data_dir / "cache/reco_mc110040.npz"); rt = TruthTable.load(cfg.data_dir / "cache/reco_mc110040_truthcols.npz")
    passed = z["passed"]
    assert int(passed.sum()) == 43643
    gr = ch.binning.digitize(z["reco_pT"][passed], z["reco_pz"][passed])
    s = ch.is_signal(rt)[passed]
    assert int((gr >= 0).sum()) == 43361 and int(((gr >= 0) & s).sum()) == 43266
    xt, yt = ch.observables(rt); gt = ch.binning.digitize(xt[passed], yt[passed])
    assert int(((gr >= 0) & (gt >= 0) & s).sum()) == 43175
    ps = ch.in_phase_space(rt)[passed]
    assert int((s & ps).sum()) == 41948 and int((s & ps & (gt >= 0)).sum()) == 41922


def test_beam_frame_sharpens_resolution():
    if not have_mc_cache():
        skip("needs the cached MC tables")
    from ndp.channels import observables as obs
    cfg = site(); ch = load_channel("minerva_me_cc_inclusive_ptpz")
    z = np.load(cfg.data_dir / "cache/reco_mc110040.npz"); rt = TruthTable.load(cfg.data_dir / "cache/reco_mc110040_truthcols.npz")
    m = z["passed"] & ch.is_signal(rt)
    def mad(a):
        med = np.median(a); return 1.4826 * np.median(np.abs(a - med))
    w_det = mad((z["reco_pT"] - obs.lep_pT(rt, "detector"))[m])
    w_beam = mad((z["reco_pT"] - obs.lep_pT(rt, "beam"))[m])
    assert w_beam < 0.08 and w_det > 0.2 and w_beam < 0.3 * w_det


def test_reference_mc_folds_to_its_own_reco():
    if not have_mc_cache():
        skip("needs the cached MC tables")
    from ndp.surrogate.base import load_surrogate
    cfg = site(); ch = load_channel("minerva_me_cc_inclusive_ptpz")
    sdir = cfg.repo_root / "surrogates/minerva_me_cc_inclusive_ptpz/binned_mc110040"
    if not (sdir / "surrogate.json").exists():
        skip("surrogate not built")
    s = load_surrogate(sdir)
    truth = TruthTable.load(cfg.data_dir / "cache/truth_mc110040.npz")
    sumw, _, _, _ = ch.truth_cells(truth)
    pred = s.fold(sumw) + s.background(truth.norm.pot)
    z = np.load(cfg.data_dir / "cache/reco_mc110040.npz")
    h, _, _ = ch.binning.histogram(z["reco_pT"][z["passed"]], z["reco_pz"][z["passed"]])
    assert np.allclose(pred, h, atol=1e-6)
