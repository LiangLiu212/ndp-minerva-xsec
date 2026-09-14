"""Registry-free per-job records.

Each grid submission is one mutable JSON record:
    <state_dir>/runs/<runtype>-YYYY-MM-DD/<stem>.gridlog

There is **no central registry**: the jobid `<runtype>-<stem>-<6hex>` decodes
itself and `find_record_for_jobid` globs for the embedded id. And there is **no
supervisor** — the work runs on the grid, so status is refreshed on demand by
polling jobsub_q (lib/monitor.py) and persisted back into the same record.

Status model (set by submit + monitor):
    pending   record written, submit not yet attempted
    submitted jobsub_submit returned a cluster id
    running   ≥1 process idle/running/transferring
    held      ≥1 process held
    done      queue drained + all outputs present
    partial   queue drained + some (< n_jobs) outputs
    failed    submit failed, or queue drained + zero outputs
    cancelled jobsub_rm called
"""
from __future__ import annotations

import json
import os
import secrets
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, Optional

from lib import config

TERMINAL_STATES = frozenset({"done", "partial", "failed", "cancelled"})
LOG_SUFFIX = ".gridlog"


def runs_root() -> Path:
    return config.state_dir() / "runs"


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _git_sha() -> Optional[str]:
    """HEAD of the *project* the state dir lives in (None outside a repo)."""
    try:
        out = subprocess.run(
            ["git", "-C", str(config.project_root()), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=False, timeout=5,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return None


def _git_dirty() -> Optional[bool]:
    """Uncommitted changes in the project worktree (None outside a repo) —
    git_sha alone cannot reproduce a dirty tree."""
    try:
        out = subprocess.run(
            ["git", "-C", str(config.project_root()), "status", "--porcelain"],
            capture_output=True, text=True, check=False, timeout=5,
        )
        if out.returncode == 0:
            return bool(out.stdout.strip())
    except Exception:
        pass
    return None


# ── jobid (self-describing; no registry) ──────────────────────────────────────

def make_jobid(runtype: str, stem: str) -> str:
    return f"{runtype}-{stem}-{secrets.token_hex(3)}"


def parse_jobid(jobid: str) -> tuple[str, str, str]:
    """Return (runtype, stem, hex6). Raises ValueError on malformed input."""
    parts = jobid.split("-")
    if len(parts) < 3:
        raise ValueError(f"malformed jobid: {jobid!r}")
    runtype = parts[0]
    hex6 = parts[-1]
    stem = "-".join(parts[1:-1])
    if not stem:
        raise ValueError(f"malformed jobid (empty stem): {jobid!r}")
    return runtype, stem, hex6


# ── paths ──────────────────────────────────────────────────────────────────────

def new_run_dir(runtype: str, when: datetime | None = None) -> Path:
    """Return (and create) `<state_dir>/runs/<runtype>-YYYY-MM-DD/`."""
    when = when or datetime.now()
    run_dir = runs_root() / f"{runtype}-{when.strftime('%Y-%m-%d')}"
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def record_path(run_dir: Path, stem: str) -> Path:
    return run_dir / f"{stem}{LOG_SUFFIX}"


# ── atomic read/write ──────────────────────────────────────────────────────────

def atomic_write_json(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    os.replace(tmp, path)


def read_record(path: Path) -> dict:
    return json.loads(Path(path).read_text())


def update_record(path: Path, **fields: Any) -> dict:
    record = read_record(path)
    record.update(fields)
    atomic_write_json(Path(path), record)
    return record


# ── schema ──────────────────────────────────────────────────────────────────────

def make_initial_record(
    *,
    jobid: str,
    runtype: str,
    stem: str,
    n_jobs: int,
    submit_user: str = "",
    command_str: str = "",
    command_file: str = "",
    submit_log_file: str = "",
    tarball_path: str = "",
    worker_script: str = "",
    pnfs_output_dir: str = "",
    local_output_dir: str = "",
    inputs: Optional[list] = None,
    outputs: Optional[dict] = None,
    extra: Optional[dict] = None,
) -> dict:
    return {
        "jobid":           jobid,
        "runtype":         runtype,
        "stem":            stem,
        "cluster_id":      "",
        "status":          "pending",
        "n_jobs":          n_jobs,
        "submitted":       utc_now_iso(),
        "finished":        None,
        "processes_done":  0,
        "processes_done_source": "",
        "outputs_pulled":  False,
        "submit_user":     submit_user,
        "command_str":     command_str,
        "command_file":    command_file,
        "submit_log_file": submit_log_file,
        "tarball_path":    tarball_path,
        "worker_script":   worker_script,
        "pnfs_output_dir": pnfs_output_dir,
        "local_output_dir": local_output_dir,
        "fetchlog_error":  "",
        "inputs":          inputs if inputs is not None else [],
        "outputs":         outputs if outputs is not None else {},
        "extra":           extra if extra is not None else {},
        "git_sha":         _git_sha(),
        "git_dirty":       _git_dirty(),
    }


# ── discovery ──────────────────────────────────────────────────────────────────

def _record_jobid(path: Path) -> Optional[str]:
    try:
        return read_record(path).get("jobid")
    except Exception:
        return None


def find_record_for_jobid(jobid: str) -> Path:
    """Resolve a jobid to its <stem>.gridlog by globbing runs/*/."""
    _, stem, _ = parse_jobid(jobid)
    matches = [p for p in runs_root().glob(f"*/{stem}{LOG_SUFFIX}")
               if _record_jobid(p) == jobid]
    if not matches:
        # stem reconstruction assumes a hyphen-free runtype; fall back to an
        # exact-jobid scan so hyphenated runtypes still resolve
        matches = [p for p in runs_root().glob(f"*/*{LOG_SUFFIX}")
                   if _record_jobid(p) == jobid]
    if not matches:
        raise FileNotFoundError(f"no record for jobid {jobid}")
    if len(matches) > 1:
        raise RuntimeError(f"multiple records for jobid {jobid}: {matches}")
    return matches[0]


def iter_record_paths() -> Iterator[Path]:
    """Yield every <stem>.gridlog under runs/, newest dir first."""
    root = runs_root()
    if not root.exists():
        return
    for d in sorted(root.iterdir(), reverse=True):
        if d.is_dir():
            yield from sorted(d.glob(f"*{LOG_SUFFIX}"))
