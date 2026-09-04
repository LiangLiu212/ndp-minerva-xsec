"""Round-trip the binning module against the paper's anc/bin_mapping.txt."""
import json
from pathlib import Path

import numpy as np
import pytest

from xsec import binning

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def mapping_rows():
    pub = json.loads((REPO_ROOT / "config" / "published.json").read_text())
    path = Path(pub["anc_dir"]) / pub["key_files"]["bin_mapping"]
    if not path.exists():
        pytest.skip(f"answer key not available at {path}")
    return binning.parse_bin_mapping(path)


def test_grid_shape():
    assert binning.N_PT == 14 and binning.N_PL == 16 and binning.N_CELLS == 224


def test_mapping_has_224_rows(mapping_rows):
    assert len(mapping_rows) == 224
    assert sorted(r["gid"] for r in mapping_rows) == list(range(224))


def test_global_id_convention(mapping_rows):
    for r in mapping_rows:
        assert binning.global_id(r["pl_bin"], r["pt_bin"]) == r["gid"]


def test_edges_match_anc_after_2dp_rounding(mapping_rows):
    # the anc file prints edges with printf("%.2f"); mimic exactly that
    # (numpy's round() differs: round(float64(0.075), 2) == 0.08, but
    # "%.2f" of the same double is "0.07")
    fmt = lambda v: float(f"{v:.2f}")
    for r in mapping_rows:
        ipl, ipt = r["pl_bin"] - 1, r["pt_bin"] - 1
        assert fmt(binning.PL_EDGES_GEV[ipl]) == r["pl_lo"]
        assert fmt(binning.PL_EDGES_GEV[ipl + 1]) == r["pl_hi"]
        assert fmt(binning.PT_EDGES_GEV[ipt]) == r["pt_lo"]
        assert fmt(binning.PT_EDGES_GEV[ipt + 1]) == r["pt_hi"]


def test_cell_ids_round_trip_bin_centers(mapping_rows):
    for r in mapping_rows:
        ipl, ipt = r["pl_bin"] - 1, r["pt_bin"] - 1
        pt_c = 0.5 * (binning.PT_EDGES_GEV[ipt] + binning.PT_EDGES_GEV[ipt + 1])
        pl_c = 0.5 * (binning.PL_EDGES_GEV[ipl] + binning.PL_EDGES_GEV[ipl + 1])
        assert binning.cell_ids([pt_c], [pl_c])[0] == r["gid"]


def test_boundary_semantics():
    # low edges inside, top edges outside ([low, high) convention)
    assert binning.cell_ids([0.0], [1.5])[0] == 0
    assert binning.cell_ids([4.5], [10.0])[0] == -1     # pT top edge
    assert binning.cell_ids([1.0], [60.0])[0] == -1     # p|| top edge
    assert binning.cell_ids([-0.1], [10.0])[0] == -1    # below pT range
    assert binning.cell_ids([1.0], [1.4])[0] == -1      # below p|| range


def test_migration_matrix_drops_outside():
    m = binning.migration_matrix([5, -1, 7], [5, 3, -1])
    assert m.sum() == 1 and m[5, 5] == 1


# ---------------------------------------------- count-conserving slots (E0)
def test_slot_scheme_dims():
    assert binning.N_PT_SLOTS == 16 and binning.N_PL_SLOTS == 18
    assert binning.N_SLOTS == 288
    assert binning.MEAS_SLOTS.shape == (224,)
    assert len(set(binning.MEAS_SLOTS.tolist())) == 224          # all distinct
    assert binning.MEAS_SLOTS.min() >= 0 and binning.MEAS_SLOTS.max() < 288


def test_slot_ids_never_drops():
    # spans every region: in-grid, pT overflow, p|| under/overflow
    pt = np.array([0.1, 2.0, 5.0, 0.1, 0.1])
    pl = np.array([5.0, 30.0, 8.0, 1.0, 70.0])
    s = binning.slot_ids(pt, pl)
    assert s.min() >= 0 and s.max() < binning.N_SLOTS
    assert (s != -1).all()


def test_meas_slots_match_cell_centers(mapping_rows):
    for r in mapping_rows:
        ipt, ipl = r["pt_bin"] - 1, r["pl_bin"] - 1
        pt_c = 0.5 * (binning.PT_EDGES_GEV[ipt] + binning.PT_EDGES_GEV[ipt + 1])
        pl_c = 0.5 * (binning.PL_EDGES_GEV[ipl] + binning.PL_EDGES_GEV[ipl + 1])
        assert binning.slot_ids([pt_c], [pl_c])[0] == binning.MEAS_SLOTS[r["gid"]]


def test_pt_underflow_impossible_pl_underflow_reachable():
    # physical p_T >= 0 -> never the pT-underflow row (slot >= N_PL_SLOTS)
    assert binning.slot_ids([0.0], [5.0])[0] >= binning.N_PL_SLOTS
    # reco p|| < 1.5 IS reachable -> pT-row but p||-underflow column (ipl=0)
    s = binning.slot_ids([0.1], [1.0])[0]
    assert s % binning.N_PL_SLOTS == 0


def test_count_conservation():
    pt = np.array([0.1, 2.0, 5.0, 0.1, 0.1, 0.3])
    pl = np.array([5.0, 30.0, 8.0, 1.0, 70.0, 6.0])
    assert binning.hist_slots(pt, pl).sum() == len(pt)            # nothing dropped
    m = binning.migration_slots(binning.slot_ids(pt, pl),
                                binning.slot_ids(pt, pl))
    assert m.sum() == len(pt)


def test_migration_colsum_is_effnum_orientation():
    # M[reco, true]; sum over axis 0 (reco) -> per-true-slot numerator
    reco = np.array([19, 20, 19])
    true = np.array([19, 19, 37])
    m = binning.migration_slots(reco, true)
    assert m[20, 19] == 1 and m[19, 19] == 1     # [reco, true]
    assert m.sum(axis=0)[19] == 2 and m.sum(axis=0)[37] == 1


def test_to_measurement_consistency_with_cell_ids():
    # for in-grid events, slot histogram projected to measurement == hist2d
    rng = np.random.default_rng(0)
    pt = rng.uniform(0.0, 4.49, 5000)
    pl = rng.uniform(1.5, 59.9, 5000)
    proj = binning.to_measurement(binning.hist_slots(pt, pl))
    flat = binning.hist2d(pt, pl).ravel()        # (14,16) row-major == GlobalID
    assert np.array_equal(proj, flat)
    assert proj.sum() == 5000                    # all in-grid here
