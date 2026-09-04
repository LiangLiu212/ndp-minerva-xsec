import numpy as np
from _helpers import toy_truth
from ndp.channels.binning import Binning
from ndp.surrogate.binned import BinnedResponse
from ndp.surrogate.parametric import SmearingSurrogate
from ndp.surrogate.base import load_surrogate

B = Binning("lep_pT", "lep_pz", (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.5), (1.5, 3.0, 4.5, 6.0, 8.0, 12.0, 20.0))


def _toy_pairs(n=60000, seed=3):
    rng = np.random.default_rng(seed)
    xt = rng.uniform(0.0, 2.5, n); yt = rng.uniform(1.5, 20.0, n)
    acc = rng.random(n) < (0.4 + 0.5 * (yt / 20.0))                 # efficiency rising with pz
    xr = xt + rng.normal(0.0, 0.08, n); yr = yt * rng.normal(1.0, 0.07, n)
    return xt, yt, xr, yr, acc


def test_binned_closure_and_persistence(tmp_path=None):
    import tempfile
    xt, yt, xr, yr, acc = _toy_pairs()
    s = BinnedResponse.fit(B, x_true_den=xt, y_true_den=yt, x_true_num=xt[acc], y_true_num=yt[acc],
                           x_reco_num=xr[acc], y_reco_num=yr[acc], pot_mc=1.0)
    true_cells, _, _ = B.histogram(xt, yt)
    pred = s.fold(true_cells) + s.background(1.0)
    reco_cells, _, _ = B.histogram(xr[acc], yr[acc])
    assert np.allclose(pred, reco_cells), "folding the training truth must reproduce the training reco exactly"
    col = s.P.sum(axis=0)
    assert np.all(col <= 1 + 1e-12)
    d = tempfile.mkdtemp()
    s.save(d); s2 = load_surrogate(d)
    assert np.allclose(s2.P, s.P) and np.allclose(s2.eff, s.eff) and s2.kind == "binned_response"
    var = s.fold_variance(true_cells)
    assert var.shape == (B.n_cells,) and np.all(var >= 0)


def test_binned_scales_linearly_and_reweights():
    xt, yt, xr, yr, acc = _toy_pairs()
    s = BinnedResponse.fit(B, x_true_den=xt, y_true_den=yt, x_true_num=xt[acc], y_true_num=yt[acc], x_reco_num=xr[acc], y_reco_num=yr[acc])
    t, _, _ = B.histogram(xt, yt)
    assert np.allclose(s.fold(2 * t), 2 * s.fold(t))
    # a pz-dependent reweight moves the folded pz projection the same direction
    w = 1 + 0.5 * (yt > 8)
    tw, _, _ = B.histogram(xt, yt, w)
    hi = B.to_grid(s.fold(tw))[:, 3:].sum(); lo = B.to_grid(s.fold(t))[:, 3:].sum()
    assert hi > lo


def test_parametric_recovers_resolution_and_folds_close_to_binned():
    xt, yt, xr, yr, acc = _toy_pairs(120000)
    s = SmearingSurrogate.fit(B, x_true_den=xt, y_true_den=yt, x_true_num=xt[acc], y_true_num=yt[acc], x_reco_num=xr[acc],
                              y_reco_num=yr[acc], mode="diff,ratio", n_samples=10)
    assert abs(np.nanmedian(s.ratio_std[:, 0]) - 0.08) < 0.015          # pT smearing width recovered
    assert abs(np.nanmedian(s.ratio_std[:, 1]) - 0.07) < 0.015          # pz fractional width recovered
    assert abs(np.nanmedian(s.ratio_mean[:, 1]) - 1.0) < 0.01
    b = BinnedResponse.fit(B, x_true_den=xt, y_true_den=yt, x_true_num=xt[acc], y_true_num=yt[acc], x_reco_num=xr[acc], y_reco_num=yr[acc])
    t, _, _ = B.histogram(xt, yt)
    f_par = s.fold_events(xt, yt, np.ones_like(xt), np.random.default_rng(5))
    f_bin = b.fold(t)
    # both predict the same total selected count within a few percent and agree cell by cell within stat
    assert abs(f_par.sum() / f_bin.sum() - 1) < 0.03
    rel = np.abs(f_par - f_bin) / np.sqrt(np.maximum(f_bin, 1))
    assert np.median(rel[f_bin > 50]) < 3.0
