"""The in-repository GENIE (pixi build) — skipped until `pixi run snapshot-genie-env` has run."""
import json
import shutil
import subprocess
from _helpers import site, skip


def _snapshot():
    p = site().repo_root / "external" / "genie_env.json"
    if not p.exists():
        skip("no external/genie_env.json (run: pixi run build-genie && pixi run snapshot-genie-env)")
    return json.loads(p.read_text())


def test_snapshot_describes_a_working_genie():
    d = _snapshot()
    env = d["env"]
    assert d["gevgen_responds"] is True
    assert env["GENIE"].endswith("external/genie/Generator")
    assert d["genie_version"].startswith("R-3_")
    assert "libEGPythia6.so" in " ".join(str(x) for x in (site().repo_root / "external/ROOTEGPythia6/install/lib").iterdir())


def test_gevgen_runs_from_the_snapshot_env():
    d = _snapshot()
    env = d["env"]
    r = subprocess.run(["gevgen", "-h"], env=env, capture_output=True, text=True, timeout=180)
    assert "gevgen" in (r.stdout + r.stderr)
    assert shutil.which("gntpc", path=env["PATH"]) is not None
