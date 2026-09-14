"""Poll jobsub_q and resolve grid-job status; persist into the .gridlog.

There is no supervisor for grid jobs (work runs remotely), so status is
recomputed on demand: `refresh_status` re-queries `jobsub_q --long`, aggregates
per-process HTCondor states, and — when the queue has drained — decides
done/partial/failed from the DONE sentinel in fetched worker logs (preferred)
or an ifdh PNFS count (fallback). Terminal states short-circuit. Ported from
genie-mcp's grid_manager, registry-free and generator-agnostic (the PNFS
count's expected suffix comes from `record['extra']['output_suffix']`).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from lib import control, records
from lib.submit_env import build_submit_env

# HTCondor JobStatus: 1=idle 2=running 3=removed 4=completed 5=held 6=transferring
_HTC_STATE_MAP = {1: "idle", 2: "running", 3: "cancelled", 4: "done", 5: "held", 6: "running"}


def parse_classad_blocks(text: str) -> list[dict]:
    """Split `jobsub_q --long` output into per-process key/value blocks."""
    blocks: list[dict] = []
    current: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            if current:
                blocks.append(current)
                current = {}
            continue
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        current[key.strip()] = val.strip().strip('"')
    if current:
        blocks.append(current)
    return blocks


def aggregate_classads(text: str) -> dict:
    """Aggregate per-process JobStatus into counts (pure; no subprocess)."""
    blocks = parse_classad_blocks(text)
    counts = {"idle": 0, "running": 0, "done": 0, "held": 0, "cancelled": 0}
    per_process = []
    for blk in blocks:
        try:
            js = int(blk.get("JobStatus", "0"))
        except ValueError:
            continue
        state = _HTC_STATE_MAP.get(js, "unknown")
        if state in counts:
            counts[state] += 1
        per_process.append({"process": blk.get("ProcId", ""), "state": state})
    return {
        "n_total": len(blocks),
        "n_idle": counts["idle"], "n_running": counts["running"],
        "n_done": counts["done"], "n_held": counts["held"],
        "n_cancelled": counts["cancelled"],
        "per_process": per_process,
        "empty": len(blocks) == 0,
    }


def query_jobsub_status(cluster_id: str, cfg: dict) -> dict:
    """Run `jobsub_q --jobid <cluster> -G <group> --long` and aggregate."""
    if not cluster_id:
        return {**aggregate_classads(""), "error": "cluster_id is empty"}
    cmd = [cfg["jobsub_q_bin"], "--jobid", cluster_id, "-G", cfg["default_group"], "--long"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=120, env=build_submit_env())
    except subprocess.TimeoutExpired:
        return {**aggregate_classads(""), "error": "jobsub_q timed out"}
    except FileNotFoundError:
        return {**aggregate_classads(""), "error": f"jobsub_q not found: {cfg['jobsub_q_bin']}"}
    out = aggregate_classads(proc.stdout or "")
    out["raw_returncode"] = proc.returncode
    return out


def count_pnfs_outputs(record: dict, cfg: dict) -> int:
    """Count <process>/ subdirs under pnfs_output_dir with ≥1 expected output.

    Generic: the expected suffix comes from `record['extra']['output_suffix']`
    (e.g. `.root`); if unset, any file counts. Best-effort; needs ifdh on PATH
    and returns 0 on any error.
    """
    pnfs_dir = record.get("pnfs_output_dir", "")
    if not pnfs_dir:
        return 0
    suffix = record.get("extra", {}).get("output_suffix")

    def _ls(path: str) -> list[str]:
        try:
            p = subprocess.run(["ifdh", "ls", path], capture_output=True,
                               text=True, timeout=120, env=build_submit_env())
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return []
        return [] if p.returncode != 0 else [l.strip() for l in (p.stdout or "").splitlines() if l.strip()]

    count = 0
    for line in _ls(pnfs_dir):
        if not line.endswith("/") or line.rstrip("/") == pnfs_dir.rstrip("/"):
            continue
        inner = _ls(line)
        if suffix:
            if any(f.endswith(suffix) for f in inner):
                count += 1
        elif any(not f.endswith("/") for f in inner):
            count += 1
    return count


def refresh_status(record_path: Path, cfg: dict) -> dict:
    """Re-poll + persist a job's status. Returns the record + a transient
    `cluster_state` summary (the latter is not persisted)."""
    record = records.read_record(record_path)
    if record["status"] in records.TERMINAL_STATES:
        return record
    cluster_id = record.get("cluster_id", "")
    if not cluster_id:
        return record  # nothing submitted to poll yet

    q = query_jobsub_status(cluster_id, cfg)
    new_status = record["status"]
    updates: dict = {}

    # Only trust the poll if jobsub_q actually succeeded: an exception path sets
    # q["error"]; a crashed jobsub_q (e.g. the OTEL KeyError) exits nonzero with
    # empty stdout, which is indistinguishable from a drained queue unless we
    # check raw_returncode. Draining on an unhealthy poll stamps a terminal
    # done/partial/failed verdict on a still-running cluster.
    query_ok = not q.get("error") and q.get("raw_returncode", 0) == 0

    if not query_ok:
        pass  # transient poll failure — keep status; q surfaced via cluster_state
    elif q["empty"]:
        n_jobs = record["n_jobs"]
        if record.get("outputs_pulled"):
            n = record.get("processes_done", 0)
            new_status = "done" if n == n_jobs else "partial"
        else:
            fl = control.fetch_log(record_path, cfg)
            if "error" not in fl and fl.get("returncode", 1) == 0:
                n = control.count_done_sentinel(fl["dest_dir"])["n_success"]
                source, fl_err = "fetchlog", ""
            else:
                n = count_pnfs_outputs(record, cfg)
                source, fl_err = "pnfs", (fl.get("error") or fl.get("stderr", ""))
            new_status = "done" if n == n_jobs else ("partial" if n > 0 else "failed")
            updates.update(processes_done=n, processes_done_source=source, fetchlog_error=fl_err)
    else:
        if q["n_held"] > 0:
            new_status = "held"
        elif q["n_running"] > 0 or q["n_idle"] > 0:
            new_status = "running"

    if new_status != record["status"]:
        updates["status"] = new_status
        if new_status in records.TERMINAL_STATES:
            updates["finished"] = records.utc_now_iso()
    if updates:
        record = records.update_record(record_path, **updates)

    out = dict(record)
    out["cluster_state"] = q
    return out


def list_jobs(cfg: dict, active_only: bool = False) -> list[dict]:
    """All records, refreshing non-terminal ones; optionally only active."""
    out: list[dict] = []
    for p in records.iter_record_paths():
        rec = records.read_record(p)
        if rec["status"] not in records.TERMINAL_STATES:
            rec = refresh_status(p, cfg)
        if active_only and rec["status"] in records.TERMINAL_STATES:
            continue
        out.append(rec)
    return out
