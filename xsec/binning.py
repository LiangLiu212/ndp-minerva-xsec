"""The paper's 2D (p_T x p_parallel) measurement grid and GlobalID mapping.

Grid: 14 muon-p_T bins x 16 p_parallel bins = 224 cells (arXiv:2106.16210).
GlobalID (0-based, anc/bin_mapping.txt): gid = (pt_bin-1)*16 + (pl_bin-1)
with 1-based bin numbers — equivalently gid = ipt*N_PL + ipl with 0-based
indices. Flat arrays in this package use exactly gid as the index.

Edge provenance: p_T edges are the tutorial's 3-decimal values
(runEventLoop.cpp:413 dansPTBins == frozen 1D chain); the anc file prints
them %.2f-rounded (0.075->0.07, 0.325->0.33, 0.475->0.47 — float64 printing,
verified 2026-06-12). p_parallel edges are exact in both.

Bin convention is ROOT-style [low, high): a value exactly on the top edge
(p_T = 4.5, p_parallel = 60) is OUTSIDE the grid (overflow).
"""
import numpy as np

PT_EDGES_GEV = np.array([0.0, 0.075, 0.15, 0.25, 0.325, 0.40, 0.475,
                         0.55, 0.70, 0.85, 1.00, 1.25, 1.50, 2.50, 4.50])
PL_EDGES_GEV = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
                         6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 40.0, 60.0])

N_PT = len(PT_EDGES_GEV) - 1   # 14
N_PL = len(PL_EDGES_GEV) - 1   # 16
N_CELLS = N_PT * N_PL          # 224


def global_id(pl_bin, pt_bin):
    """anc-convention GlobalID from 1-based (P||bin, Ptbin)."""
    return (np.asarray(pt_bin) - 1) * N_PL + (np.asarray(pl_bin) - 1)


def cell_ids(pt_gev, pl_gev):
    """Flat cell id (== GlobalID) per event, or -1 if outside the grid.

    [low, high) bins; top edges are outside.
    """
    ipt = np.digitize(np.asarray(pt_gev, dtype=np.float64), PT_EDGES_GEV) - 1
    ipl = np.digitize(np.asarray(pl_gev, dtype=np.float64), PL_EDGES_GEV) - 1
    inside = (ipt >= 0) & (ipt < N_PT) & (ipl >= 0) & (ipl < N_PL)
    return np.where(inside, ipt * N_PL + ipl, -1)


def hist2d(pt_gev, pl_gev, weights=None):
    """(N_PT, N_PL) in-grid counts; entries outside the grid are dropped."""
    h, _, _ = np.histogram2d(np.asarray(pt_gev, dtype=np.float64),
                             np.asarray(pl_gev, dtype=np.float64),
                             bins=[PT_EDGES_GEV, PL_EDGES_GEV],
                             weights=weights)
    return h


def migration_matrix(true_ids, reco_ids, weights=None):
    """(N_CELLS, N_CELLS) matrix[true_gid, reco_gid] for in-grid pairs.

    Pairs where either id is -1 (outside the grid) are NOT filled — count
    them separately for under/overflow bookkeeping.
    """
    true_ids = np.asarray(true_ids)
    reco_ids = np.asarray(reco_ids)
    both = (true_ids >= 0) & (reco_ids >= 0)
    edges = np.arange(N_CELLS + 1) - 0.5
    m, _, _ = np.histogram2d(true_ids[both], reco_ids[both],
                             bins=[edges, edges],
                             weights=None if weights is None else np.asarray(weights)[both])
    return m


# ---------------------------------------------------------------------------
# Count-conserving flat slots (for unfolding) — under/overflow included.
#
# Raw np.digitize gives, per axis: 0 = underflow, 1..N = measurement bins,
# N+1 = overflow. p_T underflow (slot ipt=0) is physically impossible
# (p_T >= 0, and PT_EDGES[0] = 0), but is kept for a uniform scheme and
# asserted empty by callers. p_parallel underflow IS reachable (reco p_par
# can be < 1.5 GeV/c — there is no reco p_z cut).
#
# Flat slot = ipt * N_PL_SLOTS + ipl, with ipt in 0..15, ipl in 0..17.
# The 224 measurement cells are the subset {ipt in 1..14, ipl in 1..16};
# MEAS_SLOTS maps each GlobalID to its flat slot. Migrations in slot space
# are oriented [reco_slot, true_slot] to match the D'Agostini unfolder.
# ---------------------------------------------------------------------------
N_PT_SLOTS = N_PT + 2   # 16: underflow + 14 bins + overflow
N_PL_SLOTS = N_PL + 2   # 18: underflow + 16 bins + overflow
N_SLOTS = N_PT_SLOTS * N_PL_SLOTS   # 288


def slot_ids(pt_gev, pl_gev):
    """Flat slot per event in [0, N_SLOTS) — never drops (count-conserving)."""
    ipt = np.digitize(np.asarray(pt_gev, dtype=np.float64), PT_EDGES_GEV)  # 0..15
    ipl = np.digitize(np.asarray(pl_gev, dtype=np.float64), PL_EDGES_GEV)  # 0..17
    return ipt * N_PL_SLOTS + ipl


def _meas_slots():
    gids = np.arange(N_CELLS)
    pt_bin = gids // N_PL + 1   # 1..14
    pl_bin = gids % N_PL + 1    # 1..16
    return pt_bin * N_PL_SLOTS + pl_bin


MEAS_SLOTS = _meas_slots()   # length 224: GlobalID -> flat slot


def to_measurement(slot_vec):
    """Extract the 224 measurement cells (GlobalID order) from a slot vector."""
    return np.asarray(slot_vec)[MEAS_SLOTS]


def hist_slots(pt_gev, pl_gev, weights=None):
    """Length-N_SLOTS histogram over flat slots; conserves all entries."""
    edges = np.arange(N_SLOTS + 1) - 0.5
    h, _ = np.histogram(slot_ids(pt_gev, pl_gev), bins=edges, weights=weights)
    return h


def migration_slots(reco_slots, true_slots, weights=None):
    """(N_SLOTS, N_SLOTS) migration M[reco_slot, true_slot], no drops.

    Oriented [reco, true] for the D'Agostini unfolder: response is the
    column-normalization over true, and M.sum(axis=0) (over reco) gives the
    efficiency numerator per true slot — free, because nothing is dropped.
    """
    edges = np.arange(N_SLOTS + 1) - 0.5
    m, _, _ = np.histogram2d(np.asarray(reco_slots), np.asarray(true_slots),
                             bins=[edges, edges], weights=weights)
    return m


def parse_bin_mapping(path):
    """Parse anc/bin_mapping.txt -> list of row dicts (224 rows)."""
    rows = []
    with open(path) as fh:
        header = fh.readline().strip().split(",")
        for line in fh:
            vals = line.strip().split(",")
            if len(vals) != len(header):
                continue
            rows.append({
                "gid": int(vals[0]), "pl_bin": int(vals[1]), "pt_bin": int(vals[2]),
                "pl_lo": float(vals[3]), "pl_hi": float(vals[4]),
                "pt_lo": float(vals[5]), "pt_hi": float(vals[6]),
            })
    return rows
