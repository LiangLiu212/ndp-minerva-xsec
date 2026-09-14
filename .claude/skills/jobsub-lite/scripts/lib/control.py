"""Cancel + fetch-log for grid jobs, plus the DONE-sentinel completion count.

`cancel`  → jobsub_rm; `fetch_log` → jobsub_fetchlog --unzipdir <stem>.fetched/;
`count_done_sentinel` counts fetched worker `*.out` logs whose body contains a
standalone `DONE` line (the authoritative "this process finished" signal emitted
by the worker scripts — more reliable than ifdh PNFS counting). All jobsub calls
run under the scrubbed env.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

from lib import records
from lib.submit_env import build_submit_env


def _jobsub_rm_bin(cfg: dict) -> str:
    return cfg.get("jobsub_rm_bin") or str(Path(cfg["jobsub_q_bin"]).parent / "jobsub_rm")


def cancel(record_path: Path, cfg: dict) -> dict:
    """jobsub_rm the job's cluster and mark the record cancelled — only when
    the rm succeeded. A failed rm keeps the prior status: cancelled is
    terminal, so stamping it would freeze refresh_status on a cluster that is
    still running."""
    record = records.read_record(record_path)
    cluster_id = record.get("cluster_id", "")
    if not cluster_id:
        rec = records.update_record(record_path, status="cancelled",
                                    finished=records.utc_now_iso())
        rec["_cancel_ok"] = True
        rec["_cancel_detail"] = "no cluster id; record marked cancelled"
        return rec

    cmd = [_jobsub_rm_bin(cfg), "--jobid", cluster_id, "-G", cfg["default_group"]]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=120, env=build_submit_env())
        ok = proc.returncode == 0
        detail = (proc.stderr or proc.stdout or "")[-500:]
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        ok, detail = False, f"{cmd[0]} failed: {e}"

    if ok:
        rec = records.update_record(record_path, status="cancelled",
                                    finished=records.utc_now_iso())
    else:
        rec = records.update_record(
            record_path,
            cancel_error=f"jobsub_rm failed; status kept: {detail}"[:500])
    rec["_cancel_ok"] = ok
    rec["_cancel_detail"] = detail
    rec["_cancel_command"] = " ".join(cmd)
    return rec


def fetch_log(record_path: Path, cfg: dict, dest_dir: Optional[str] = None) -> dict:
    """jobsub_fetchlog the job's worker logs into <stem>.fetched/ (or dest_dir)."""
    record = records.read_record(record_path)
    cluster_id = record.get("cluster_id", "")
    if not cluster_id:
        return {"error": "job has no cluster id", "returncode": 1}

    run_dir = Path(record_path).parent
    dest = Path(dest_dir) if dest_dir else run_dir / f"{record['stem']}.fetched"
    dest.mkdir(parents=True, exist_ok=True)

    cmd = [cfg["jobsub_fetchlog_bin"], "--jobid", cluster_id,
           "-G", cfg["default_group"], "--unzipdir", str(dest)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=int(cfg.get("fetchlog_timeout_s", 600)),
                              env=build_submit_env())
        rc, stderr = proc.returncode, (proc.stderr or "")[-500:]
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        return {"error": f"{cmd[0]} failed: {e}", "returncode": -1, "dest_dir": str(dest)}

    files = sorted(str(p) for p in dest.rglob("*") if p.is_file())
    return {
        "dest_dir": str(dest), "n_files": len(files), "files": files,
        "command_str": " ".join(cmd), "returncode": rc, "stderr": stderr,
    }


def count_done_sentinel(dest_dir: str | Path) -> dict:
    """Count fetched `*.out` worker logs whose body has a standalone `DONE`
    line (rglob: fetchlog layouts may nest per-process subdirs)."""
    dest = Path(dest_dir)
    n_success = n_inspected = 0
    for out in dest.rglob("*.out"):
        n_inspected += 1
        try:
            text = out.read_text(errors="replace")
        except OSError:
            continue
        if any(line.strip() == "DONE" for line in text.splitlines()):
            n_success += 1
    return {"n_success": n_success, "n_inspected": n_inspected}
