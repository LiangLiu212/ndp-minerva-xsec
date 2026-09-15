"""The chunked surrogate builder: summed count_pairs over chunks == fit on the whole sample, bit for bit;
persistence of the feed-in and per-category background arrays; equality with the one-shot builder on
the legacy single-file caches when they exist."""
import tempfile
from pathlib import Path

import numpy as np

from ndp.channels.binning import Binning
from ndp.surrogate.base import load_surrogate
from ndp.surrogate.binned import BinnedResponse, add_counts, count_pairs
from _helpers import legacy_channel, have_mc_cache, have_data_file, site, skip


def _sample(n, rng):
    xt = rng.uniform(0.0, 12.0, n); yt = rng.uniform(-1.5, 1.5, n)
    xr = xt + rng.normal(0, 0.4, n); yr = yt + rng.normal(0, 0.1, n)
    cat = rng.choice(["bkg 1 pi+-", "bkg 1 pi0", "bkg multi-pi", "bkg other (no pion)"], size=n)
    return xt, yt, xr, yr, cat


def test_chunked_counts_equal_whole_fit():
    rng = np.random.default_rng(5)
    b = Binning("x", "y", (0.0, 2.0, 5.0, 10.0), (-1.0, 0.0, 1.0))
    cats = ("bkg 1 pi+-", "bkg 1 pi0", "bkg multi-pi", "bkg other (no pion)")
    den = [_sample(3000, rng) for _ in range(3)]          # three chunks of "truth signal"
    num = [tuple(a[:900] for a in d) for d in den]         # a subset is selected
    bkg = [_sample(700, rng) for _ in range(3)]
    acc = None
    for d, nm, bk in zip(den, num, bkg):
        acc = add_counts(acc, count_pairs(b, x_true_den=d[0], y_true_den=d[1], x_true_num=nm[0], y_true_num=nm[1],
                                          x_reco_num=nm[2], y_reco_num=nm[3], x_reco_bkg=bk[2], y_reco_bkg=bk[3],
                                          bkg_category=bk[4], category_names=cats))
    whole = BinnedResponse.fit(b, x_true_den=np.concatenate([d[0] for d in den]), y_true_den=np.concatenate([d[1] for d in den]),
                               x_true_num=np.concatenate([n[0] for n in num]), y_true_num=np.concatenate([n[1] for n in num]),
                               x_reco_num=np.concatenate([n[2] for n in num]), y_reco_num=np.concatenate([n[3] for n in num]),
                               x_reco_bkg=np.concatenate([k[2] for k in bkg]), y_reco_bkg=np.concatenate([k[3] for k in bkg]), pot_mc=2.0e20)
    chunked = BinnedResponse.from_counts(b, acc, pot_mc=2.0e20)
    for k in ("eff", "P", "migration_counts", "den_counts", "num_counts", "bkg_per_pot"):
        assert np.array_equal(getattr(whole, k), getattr(chunked, k)), k
    # the category split sums to the explicit background; feed-in + explicit == the stored total
    by_cat = sum(chunked.bkg_by_category_counts[c] for c in cats)
    assert np.array_equal(by_cat, acc["bkg_explicit"])
    assert np.allclose(chunked.feedin_counts + acc["bkg_explicit"], chunked.bkg_per_pot * 2.0e20)
    assert np.allclose(chunked.fold_eff_only(chunked.den_counts), chunked.num_counts)
    # closure: fold(den) + bkg == all selected reco events inside the grid
    sel_x = np.concatenate([n[2] for n in num] + [k[2] for k in bkg]); sel_y = np.concatenate([n[3] for n in num] + [k[3] for k in bkg])
    h, _, _ = b.histogram(sel_x, sel_y)
    assert np.allclose(chunked.fold(chunked.den_counts) + chunked.background(2.0e20), h)
    # persistence keeps the new arrays
    d = Path(tempfile.mkdtemp()) / "sur"
    chunked.save(d)
    back = load_surrogate(d)
    assert np.array_equal(back.feedin_counts, chunked.feedin_counts)
    assert set(back.bkg_by_category_counts) == set(cats) and np.array_equal(back.bkg_by_category_counts[cats[1]], chunked.bkg_by_category_counts[cats[1]])
    assert back.background_by_category(1.0e20)[cats[0]].sum() == chunked.bkg_by_category_counts[cats[0]].sum() / 2.0


def test_chunked_builder_matches_one_shot_on_legacy_caches():
    if not (have_mc_cache() and have_data_file()):
        skip("needs the cached single-file MC tables")
    from ndp.channels import load_channel, load_measurement
    from ndp.surrogate.build import build_surrogates
    from ndp.surrogate.chunked import build_surrogates_chunked
    cfg = site(); ch = legacy_channel("minerva_me_ccqelike_1mu1p"); m = load_measurement(ch, "muon_p")
    root = Path(tempfile.mkdtemp())
    one = build_surrogates(ch, m, cfg, kinds=("binned",), out_root=root / "one", log=lambda *a: None)[0]
    chk = build_surrogates_chunked(ch, [m], cfg, out_root=root / "chunk", log=lambda *a: None)["muon_p"]
    a, b = load_surrogate(one["path"]), load_surrogate(chk["path"])
    for k in ("eff", "P", "migration_counts", "den_counts", "num_counts", "bkg_per_pot"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k
    assert chk["closure"]["exact"] and one["closure"]["exact"]
