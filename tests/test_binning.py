import numpy as np
from _helpers import site, have_minerva_repo, skip
from ndp.channels.binning import Binning

PT = (0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5, 4.5)
PZ = (1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 40.0, 60.0)


def test_cell_roundtrip_both_formulas():
    for f in ("ix*n_y + iy", "iy*n_x + ix", "ipt*n_pz + ipz", "ipz*n_pt + ipt"):
        b = Binning("lep_pT", "lep_pz", PT, PZ, f)
        ix, iy = b.unravel()
        assert np.array_equal(b.cell(ix, iy), np.arange(b.n_cells))
        g = b.to_grid(np.arange(b.n_cells))
        assert np.array_equal(b.from_grid(g), np.arange(b.n_cells))


def test_digitize_edges_and_outside():
    b = Binning("x", "y", PT, PZ)
    g = b.digitize(np.array([0.0, 0.075, 4.5, 2.0, -0.1]), np.array([1.5, 1.5, 3.0, 60.0, 3.0]))
    assert g[0] == b.cell(0, 0) and g[1] == b.cell(1, 0)
    assert g[2] == -1 and g[3] == -1 and g[4] == -1     # upper edges exclusive, negatives outside


def test_histogram_and_projection_conserve():
    rng = np.random.default_rng(1)
    b = Binning("x", "y", PT, PZ)
    x, y = rng.uniform(0, 5, 5000), rng.uniform(1, 70, 5000)
    h, h2, n_out = b.histogram(x, y)
    inside = (x < 4.5) & (y >= 1.5) & (y < 60)
    assert h.sum() == inside.sum() and n_out == (~inside).sum()
    assert abs(b.project(h, "x", per_width=False).sum() - h.sum()) < 1e-9
    assert np.allclose(b.areas().sum(), 4.5 * 58.5)


def test_matches_published_bin_mapping():
    if not have_minerva_repo():
        skip("MINERvA repo not configured")
    path = site().minerva_repo / "papers/minerva/2106.16210/anc/bin_mapping.txt"
    if not path.exists():
        skip("bin_mapping.txt not present")
    b = Binning("lep_pT", "lep_pz", PT, PZ, "ipt*n_pz + ipz")
    import csv
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 224
    for r in rows:
        g = int(r["GlobalID"]); ipl = int(r["P||bin"]) - 1; ipt = int(r["Ptbin"]) - 1
        assert b.cell(ipt, ipl) == g
