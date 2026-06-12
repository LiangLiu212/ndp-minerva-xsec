"""Offline consistency tests for the Stage-1 config artifacts.

No network: validates the dataset specs, branch catalog, and run-period table
against each other and against the certified constants they mirror.
"""
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
CFG = REPO_ROOT / "config"

# The certified reco-selection branch set (tools/cc_inclusive_selector.py,
# exploration repo). The catalog's reco_selection role must match exactly.
CERTIFIED_RECO_BRANCHES = {
    "vtx",
    "muon_thetaX",
    "muon_thetaY",
    "isMinosMatchTrack",
    "phys_n_dead_discr_pair_upstream_prim_track_proj",
    "MasterAnaDev_minos_trk_qp",
}

OPEN_DATA_PLAYLISTS = 12  # 1A-1G, 1L-1P (1H is transitional, not released)


@pytest.fixture(scope="module")
def catalog():
    return json.loads((CFG / "branches.json").read_text())


@pytest.fixture(scope="module")
def spec():
    return json.loads((CFG / "datasets" / "me1A_single_pair.json").read_text())


@pytest.fixture(scope="module")
def periods():
    return json.loads((CFG / "run_periods.json").read_text())


def test_all_configs_parse():
    for p in [*CFG.glob("*.json"), *(CFG / "datasets").glob("*.json")]:
        json.loads(p.read_text())


def test_catalog_reco_selection_matches_certified(catalog):
    got = {b["name"] for b in catalog["branches"] if "reco_selection" in b["roles"]}
    assert got == CERTIFIED_RECO_BRANCHES


def test_catalog_roles_are_declared(catalog):
    declared = set(catalog["_meta"]["roles"])
    used = {r for b in catalog["branches"] for r in b["roles"]}
    assert used <= declared, f"undeclared roles: {used - declared}"


def test_catalog_required_role_coverage(catalog):
    used = {r for b in catalog["branches"] for r in b["roles"]}
    for needed in ("measurement", "measurement_true", "reco_selection",
                   "signal_definition", "phase_space", "pot"):
        assert needed in used, f"no branch carries role {needed}"


def test_catalog_mc_only_branches_marked_absent_in_data(catalog):
    for b in catalog["branches"]:
        if b["name"].startswith("mc_") or b.get("family"):
            assert b["in_data"] is False, f"{b['name']} must be absent in data"


def test_spec_urls_come_from_official_lists(spec):
    lists_dir = CFG / "playlists"
    official = set()
    for lst in lists_dir.glob("*.txt"):
        official |= {ln.strip() for ln in lst.read_text().splitlines() if ln.strip()}
    for entry in spec["files"]:
        assert entry["url"] in official, f"{entry['url']} not in official playlist lists"


def test_spec_expected_fingerprints_present(spec):
    for entry in spec["files"]:
        exp = entry["expected"]
        for key in ("trees", "meta_entries", "pot_used", "pot_total", "reco_entries"):
            assert key in exp, f"{entry['role']} spec missing expected.{key}"
        if entry["role"] == "mc":
            assert "truth_entries" in exp


def test_run_periods_open_data_count(periods):
    pls = periods["playlists"]
    assert sum(1 for p in pls.values() if p["in_open_data"]) == OPEN_DATA_PLAYLISTS
    assert pls["minervame1H"]["in_open_data"] is False


def test_run_periods_me1A_window(periods):
    p = periods["playlists"]["minervame1A"]
    assert p["run_first"] == 6038 and p["run_last"] == 10066
    assert p["t_start"].startswith("2013-09-12") and p["t_end"].startswith("2014-01-14")


def test_run_periods_pot_sum_matches_getdata_page(periods):
    total = sum(p["data_pot_e20"] or 0 for p in periods["playlists"].values())
    assert abs(total - 11.12) < 0.005, f"sum {total} != 11.12e20 (getdata page)"


def test_golden_runs_inside_their_playlist_window(spec, periods):
    p = periods["playlists"][spec["playlist"]]
    for entry in spec["files"]:
        if entry["role"] == "data":
            assert p["run_first"] <= entry["run"] <= p["run_last"]
