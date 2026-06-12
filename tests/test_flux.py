"""Tests for xsec.flux (flux CV weight, MnvTune v1 weighter #1).

Need the unpacked FluxAndReweightFiles under data/flux/ — skipped otherwise
(local files; no network). The headline gate: the nu-e-constrained flux
integral must reproduce the paper's 6.32e-8 nu/cm^2/POT within 5%.
"""
import numpy as np
import pytest

from xsec import flux

pytestmark = pytest.mark.skipif(
    not flux.FLUX_DIR.exists() or not flux.CONSTRAINT_FILE.exists(),
    reason="FluxAndReweightFiles not unpacked under data/flux/")

PAPER_PHI_CM2 = 6.32e-8  # nu/cm^2/POT, 0-100 GeV (arXiv:2106.16210)


@pytest.fixture(scope="module")
def fx():
    return flux.FluxCV("minervame1A")


def test_playlist_group_mapping():
    assert flux.FLUX_GROUP["minervame1A"] == "minervame1D"
    assert flux.FLUX_GROUP["minervame1G"] == "minervame1M"
    assert flux.FLUX_GROUP["minervame1P"] == "minervame1N"
    assert len(flux.FLUX_GROUP) == 12


def test_constraint_weights_normalized():
    w = flux._load_constraint_weights("Flux")
    assert len(w) == 1000
    assert abs(w.sum() - 1000.0) < 1e-6  # file ships normalized to N


def test_universe_count(fx):
    assert fx.n_universes == 1000


def test_phi_integral_gate(fx):
    phi_cm2 = fx.integral(0.0, 100.0) * 1e-4
    assert abs(phi_cm2 / PAPER_PHI_CM2 - 1) < 0.05, f"phi={phi_cm2:.4e}"


def test_constraint_pulls_flux_down(fx):
    con = fx.integral(which="constrained")
    rew = fx.integral(which="reweighted")
    assert con < rew  # nu-e data pulls the PPFX flux down
    assert 0.85 < con / rew < 0.95  # observed 0.901
    real = fx.correction[1:-1]
    assert real.min() > 0.5 and real.max() < 1.2


def test_cv_weight_matches_root_interpolate(fx):
    ROOT = pytest.importorskip("ROOT")
    nb = len(fx.centers)
    h_con = ROOT.TH1D("hcon", "", nb, fx.edges)
    h_gen = ROOT.TH1D("hgen", "", nb, fx.edges)
    for i in range(nb):
        h_con.SetBinContent(i + 1, fx._con_real[i])
        h_gen.SetBinContent(i + 1, fx._gen_real[i])
    rng = np.random.default_rng(42)
    enus = np.concatenate([rng.uniform(0.1, 74.9, 200),
                           rng.uniform(75.1, 99.9, 50)])
    ours = fx.cv_weight(enus)
    for e, w in zip(enus, ours):
        if e > flux.INTERPOLATE_MAX_GEV:
            num = h_con.GetBinContent(h_con.FindBin(e))
            den = h_gen.GetBinContent(h_gen.FindBin(e))
        else:
            num, den = h_con.Interpolate(e), h_gen.Interpolate(e)
        expect = num / den if num != 0 and den != 0 else 1.0
        assert abs(w - expect) <= 1e-12 * max(1.0, abs(expect)), f"Enu={e}"


def test_cv_weight_zero_guard(fx):
    # far outside any populated bin the clamped interpolation still returns
    # finite content, but a synthetic zero must give weight 1
    w = fx.cv_weight(np.array([5.0]))
    assert np.isfinite(w).all()
    saved = fx._gen_real.copy()
    try:
        fx._gen_real[:] = 0.0
        assert fx.cv_weight(np.array([5.0]))[0] == 1.0
    finally:
        fx._gen_real[:] = saved


def test_cv_weight_envelope(fx):
    # Observed envelope (regression band): the weight pulls the focusing-peak
    # region DOWN (~0.75-0.88 below ~10 GeV) and the falling tail UP (to ~1.64
    # near 45 GeV) — the shape that explains the unweighted data/MC ratios
    # (low-p|| ~0.85, 40-60 GeV p|| ~1.9). See docs/cv_reweight.md §1.
    enus = np.linspace(1.5, 60.0, 500)
    w = fx.cv_weight(enus)
    assert np.isfinite(w).all()
    assert 0.70 < w.min() < 0.80, w.min()
    assert 1.55 < w.max() < 1.70, w.max()
    # anchor points of the measured curve (non-monotonic: bump at ~15 GeV,
    # dip at ~25 GeV, peak ~1.61 near 50 GeV)
    w_at = lambda e: fx.cv_weight(np.array([e]))[0]
    assert 0.78 < w_at(5.0) < 0.85
    assert w_at(50.0) > 1.4
