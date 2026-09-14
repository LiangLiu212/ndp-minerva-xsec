"""Per-project state dir (`.jobsub/`) + config.json loading and generation.

All mutable state — config.json, catalog.json, tarballs/, runs/ — lives in a
`.jobsub/` directory inside the *project*, never inside the skill. Resolution
order for the state dir: `$JOBSUB_STATE_DIR` → nearest ancestor of CWD (CWD
included) containing `.jobsub/` → error telling the user to run `init`.

Config path precedence: explicit `path=` arg → `$JOBSUB_CONFIG` env →
`<state_dir>/config.json`.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

STATE_DIR_NAME = ".jobsub"

_REQUIRED_KEYS = (
    "jobsub_bin",
    "jobsub_q_bin",
    "jobsub_fetchlog_bin",
    "default_group",
)

_JOBSUB_BIN_DIR = Path("/opt/jobsub_lite/bin")
_JOBSUB_BINS = ("jobsub_submit", "jobsub_q", "jobsub_rm", "jobsub_fetchlog")


# ── state dir ──────────────────────────────────────────────────────────────────

def find_state_dir(start: Path | None = None) -> Path | None:
    """Nearest `<ancestor>/.jobsub` of `start` (default CWD), or None."""
    env = os.environ.get("JOBSUB_STATE_DIR")
    if env:
        p = Path(env)
        if not p.is_dir():
            # A typo'd override would otherwise silently grow a fresh state
            # tree at the bogus path.
            raise FileNotFoundError(
                f"$JOBSUB_STATE_DIR points at a nonexistent dir: {p} — "
                f"unset it, or create the state dir with `jobsub.py init`")
        return p
    d = (start or Path.cwd()).resolve()
    for candidate in (d, *d.parents):
        p = candidate / STATE_DIR_NAME
        if p.is_dir():
            return p
    return None


def state_dir(start: Path | None = None) -> Path:
    p = find_state_dir(start)
    if p is None:
        raise FileNotFoundError(
            f"no {STATE_DIR_NAME}/ found in {Path.cwd()} or any parent "
            f"(and $JOBSUB_STATE_DIR unset) — run `jobsub.py init` from the "
            f"project root first"
        )
    return p


def project_root(start: Path | None = None) -> Path:
    return state_dir(start).parent


def ensure_state_dir(base: Path) -> Path:
    """Create `<base>/.jobsub/` (+ a self-ignoring .gitignore) and return it."""
    p = Path(base).resolve() / STATE_DIR_NAME
    p.mkdir(parents=True, exist_ok=True)
    gitignore = p / ".gitignore"
    if not gitignore.exists():
        gitignore.write_text("# jobsub-lite skill state — machine-local, never commit\n*\n")
    return p


# ── config ─────────────────────────────────────────────────────────────────────

def config_path(path: str | Path | None = None) -> Path:
    return Path(path or os.environ.get("JOBSUB_CONFIG") or state_dir() / "config.json")


def load_config(path: str | Path | None = None) -> dict:
    cfg_path = config_path(path)
    if not cfg_path.exists():
        raise FileNotFoundError(
            f"jobsub config not found: {cfg_path} "
            f"(run `jobsub.py init`, or set $JOBSUB_CONFIG)"
        )
    cfg = json.loads(cfg_path.read_text())

    missing = [k for k in _REQUIRED_KEYS if not cfg.get(k)]
    if missing:
        raise KeyError(f"{cfg_path} is missing required key(s): {', '.join(missing)}")

    # jobsub_rm is derived from jobsub_q's dir if absent, to avoid a stale path.
    cfg.setdefault("jobsub_rm_bin", str(Path(cfg["jobsub_q_bin"]).parent / "jobsub_rm"))
    cfg["config_path"] = str(cfg_path)
    return cfg


def probe_jobsub_bins() -> dict:
    """Locate the jobsub_lite binaries: /opt/jobsub_lite/bin first, then PATH.
    Values fall back to the bare name (with a warning from the caller) so the
    generated config is editable rather than unwritable."""
    out = {}
    for name in _JOBSUB_BINS:
        cand = _JOBSUB_BIN_DIR / name
        if cand.exists():
            out[f"{name}_bin" if name != "jobsub_submit" else "jobsub_bin"] = str(cand)
            continue
        found = shutil.which(name)
        key = "jobsub_bin" if name == "jobsub_submit" else f"{name}_bin"
        out[key] = found or name
    return out


def default_config(group: str) -> dict:
    cfg = probe_jobsub_bins()
    cfg.update({
        "default_group": group,
        "default_role": "Analysis",
        "default_disk": "20GB",
        "pnfs_scratch_base": f"/pnfs/{group}/scratch/users",
        "submit_timeout_s": 600,
        "fetchlog_timeout_s": 600,
    })
    return cfg


def resolve_group(explicit: str | None = None) -> str | None:
    """--group arg → $JOBSUB_GROUP → $EXPERIMENT → None (caller must ask)."""
    return explicit or os.environ.get("JOBSUB_GROUP") or os.environ.get("EXPERIMENT") or None
