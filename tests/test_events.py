import numpy as np
from _helpers import toy_truth
from ndp.events import TruthTable


def test_roundtrip_and_select(tmp_path=None):
    import tempfile, pathlib
    t = toy_truth(500, with_fs=True)
    d = pathlib.Path(tempfile.mkdtemp())
    p = t.save(d / "t.npz")
    t2 = TruthTable.load(p)
    assert t2.n == 500 and t2.has_fs
    assert np.array_equal(t2["fs_pdg"], t["fs_pdg"])
    m = t["E_nu"] > np.median(t["E_nu"])
    s = t.select(m)
    assert s.n == int(m.sum())
    # jagged selection keeps per-event particle counts
    c_full = np.diff(t["fs_offsets"])[m]
    assert np.array_equal(np.diff(s["fs_offsets"]), c_full)
    assert s["fs_offsets"][-1] == len(s["fs_pdg"])


def test_fs_sum_matches_loop():
    t = toy_truth(300, with_fs=True)
    tot = t.fs_sum(t["fs_E"], t["fs_pdg"] == 2212)
    off = t["fs_offsets"]
    for i in range(0, 300, 37):
        sl = slice(off[i], off[i + 1])
        expect = t["fs_E"][sl][t["fs_pdg"][sl] == 2212].sum()
        assert abs(tot[i] - expect) < 1e-12


def test_with_weights_and_concat():
    t = toy_truth(100)
    t2 = t.with_weights(np.full(100, 2.0), note="x2")
    assert t2["weight"].sum() == 200 and t["weight"].sum() == 100
    c = TruthTable.concatenate([toy_truth(50, with_fs=True), toy_truth(70, seed=1, with_fs=True)])
    assert c.n == 120 and c["fs_offsets"][-1] == len(c["fs_pdg"])
