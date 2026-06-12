"""Input tests: official minervame1A playlist file lists (Data + StandardMC).

Downloads the official per-playlist text files (one xrootd URL per line) published
on https://minerva.fnal.gov/getdata/ into config/playlists/ — only if not already
present (delete a file to force a re-download) — then validates their contents
against the published ME1A run ranges (https://minerva.fnal.gov/data-run-periods/).

These lists are the source of truth for which AnaTuple files exist in the open
data release; the Stage-1 dataset specs derive per-file xrootd URLs from them.
"""
import re
import urllib.error
import urllib.request
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PLAYLIST_DIR = REPO_ROOT / "config" / "playlists"

BASE_URL = "https://minerva.fnal.gov/wp-content/uploads/2026/03"
FILE_LISTS = {
    "data": "MediumEnergy_FHC_Data_Playlist1A.txt",
    "mc": "MediumEnergy_FHC_StandardMC_Playlist1A.txt",
}

# Expected line format, anchored to the documented open-data layout.
LINE_RE = {
    "data": re.compile(
        r"^root://fndcadoor\.fnal\.gov:1095/pnfs/fnal\.gov/usr/minerva/persistent/"
        r"OpenData/MediumEnergy_FHC/Data/Playlist1A/"
        r"MasterAnaDev_data_AnaTuple_run(\d{8})_Playlist\.root$"
    ),
    "mc": re.compile(
        r"^root://fndcadoor\.fnal\.gov:1095/pnfs/fnal\.gov/usr/minerva/persistent/"
        r"OpenData/MediumEnergy_FHC/MC/StandardMC/Playlist1A/"
        r"MasterAnaDev_mc_AnaTuple_run(\d{8})_Playlist\.root$"
    ),
}

# minervame1A data run range per minerva.fnal.gov/data-run-periods/ (6038/31-10066/23).
DATA_RUN_RANGE = (6038, 10066)
# StandardMC Playlist1A uses the 110xxx run block (110000...).
MC_RUN_RANGE = (110000, 119999)

# Pin the list sizes observed at first ingestion (2026-06-12). A failure here means
# the upstream release changed — re-inspect before accepting.
EXPECTED_N_FILES = {"data": 253, "mc": 41}

# Runs of the certified single-run pair (golden counts 844/43643 were derived on
# these), which must be present in the official lists.
GOLDEN_RUNS = {"data": 10066, "mc": 110040}


def _fetch(kind: str) -> Path:
    """Return the local playlist file, downloading it first if absent."""
    PLAYLIST_DIR.mkdir(parents=True, exist_ok=True)
    dest = PLAYLIST_DIR / FILE_LISTS[kind]
    if dest.exists():
        return dest
    url = f"{BASE_URL}/{FILE_LISTS[kind]}"
    req = urllib.request.Request(url, headers={"User-Agent": "ndp-minerva-xsec-tests"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            body = resp.read()
    except urllib.error.HTTPError as err:
        pytest.fail(f"HTTP {err.code} fetching {url} — URL scheme changed upstream?")
    except (urllib.error.URLError, TimeoutError) as err:
        pytest.skip(f"network unavailable for {url}: {err}")
    dest.write_bytes(body)
    return dest


@pytest.fixture(scope="module", params=["data", "mc"])
def playlist(request):
    kind = request.param
    path = _fetch(kind)
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    return kind, path, lines


def test_list_present_and_nonempty(playlist):
    kind, path, lines = playlist
    assert path.is_file(), f"{kind} playlist file missing at {path}"
    assert lines, f"{kind} playlist file {path} is empty"


def test_expected_file_count(playlist):
    kind, path, lines = playlist
    assert len(lines) == EXPECTED_N_FILES[kind], (
        f"{kind} list has {len(lines)} entries, expected {EXPECTED_N_FILES[kind]} "
        f"(pinned 2026-06-12) — upstream release changed, re-inspect {path}"
    )


def test_every_line_is_wellformed_xrootd_url(playlist):
    kind, _, lines = playlist
    bad = [ln for ln in lines if not LINE_RE[kind].match(ln)]
    assert not bad, f"{kind} list has {len(bad)} malformed lines, first: {bad[0]!r}"


def test_no_duplicate_entries(playlist):
    kind, _, lines = playlist
    assert len(set(lines)) == len(lines), f"{kind} list contains duplicate entries"


def test_run_numbers_within_published_range(playlist):
    kind, _, lines = playlist
    lo, hi = DATA_RUN_RANGE if kind == "data" else MC_RUN_RANGE
    runs = sorted(int(LINE_RE[kind].match(ln).group(1)) for ln in lines)
    assert len(set(runs)) == len(runs), f"{kind} list repeats a run number"
    assert runs[0] >= lo and runs[-1] <= hi, (
        f"{kind} runs span [{runs[0]}, {runs[-1]}], outside published [{lo}, {hi}]"
    )


def test_golden_run_present(playlist):
    kind, _, lines = playlist
    golden = GOLDEN_RUNS[kind]
    runs = {int(LINE_RE[kind].match(ln).group(1)) for ln in lines}
    assert golden in runs, (
        f"certified {kind} run {golden} not in the official {kind} list — "
        f"golden parity counts would have no anchor"
    )
