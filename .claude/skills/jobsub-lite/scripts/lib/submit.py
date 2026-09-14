"""Run a prepared jobsub_submit command and persist a per-job record.

Generic: the caller hands a **complete** `submit_cmd` (argv starting with the
jobsub_submit binary and ending with `file://<worker>` + the worker's own args).
`submit()` appends `--no_submit` for dry runs, runs it under the scrubbed env
(lib/submit_env), captures the combined log, parses the cluster id, and writes
three sibling files in the run dir:

    <stem>.command.json   the resolved argv + metadata (written first)
    <stem>.submit.log     combined jobsub_submit stdout+stderr
    <stem>.gridlog        the mutable record (pending → submitted/failed)
"""
from __future__ import annotations

import json
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

from lib import records
from lib.submit_env import build_submit_env

# <cluster>.<proc>@<host>.fnal.gov — take the LAST match in jobsub_submit output.
_RE_CLUSTER_ID = re.compile(r"\b(\d+\.\d+@\S+\.fnal\.gov)\b")


def parse_cluster_id(combined_output: str) -> str:
    matches = _RE_CLUSTER_ID.findall(combined_output)
    if not matches:
        raise ValueError(
            "no cluster id (<n>.<n>@host.fnal.gov) found in jobsub_submit output"
        )
    return matches[-1]


def submit(
    *,
    runtype: str,
    stem: str,
    submit_cmd: list[str],
    n_jobs: int,
    submit_user: str = "",
    worker_script: str = "",
    tarball_path: str = "",
    pnfs_output_dir: str = "",
    local_output_dir: str = "",
    inputs: Optional[list] = None,
    outputs: Optional[dict] = None,
    extra: Optional[dict] = None,
    dry_run: bool = False,
    exe_index: Optional[int] = None,
    timeout: int = 600,
    when: Optional[datetime] = None,
) -> tuple[dict, Path]:
    """Submit `submit_cmd`; return (record, gridlog_path).

    `exe_index` is the position of the `file://` executable inside
    `submit_cmd` (required for dry runs — `--no_submit` is inserted there).
    Raises FileExistsError when the target gridlog already holds a
    non-pending record: same runtype + stem + day would silently clobber the
    only record of a possibly-running cluster.
    """
    when = when or datetime.now()
    run_dir = records.new_run_dir(runtype, when=when)
    gridlog = records.record_path(run_dir, stem)
    if gridlog.exists():
        try:
            prev = records.read_record(gridlog)
        except Exception:
            prev = {}
        # Overwriting a pending record is the normal dry-run → real-submit
        # flow; anything else is a live/finished record being clobbered.
        if prev.get("status") != "pending":
            raise FileExistsError(
                f"record already exists for runtype='{runtype}' stem='{stem}' today "
                f"(jobid {prev.get('jobid', '?')}, status {prev.get('status', '?')}): "
                f"{gridlog} — pick a different --stem")
    command_file = run_dir / f"{stem}.command.json"
    submit_log = run_dir / f"{stem}.submit.log"

    cmd = list(submit_cmd)
    if dry_run:
        # --no_submit must precede the executable: jobsub_submit stops parsing
        # its own options at the file:// argument, so anything appended after it
        # becomes a worker argument and the job is submitted for real. The
        # executable's position comes from the caller — a first-file://
        # heuristic breaks when a `-f file:///...` input precedes it.
        if exe_index is None or not (0 <= exe_index < len(cmd)) \
                or not cmd[exe_index].startswith("file://"):
            raise ValueError(
                f"dry_run=True needs exe_index pointing at the file:// "
                f"executable in submit_cmd (got exe_index={exe_index!r})")
        cmd.insert(exe_index, "--no_submit")
    jobid = records.make_jobid(runtype, stem)

    record = records.make_initial_record(
        jobid=jobid, runtype=runtype, stem=stem, n_jobs=n_jobs,
        submit_user=submit_user, command_str=" ".join(cmd),
        command_file=str(command_file), submit_log_file=str(submit_log),
        tarball_path=tarball_path, worker_script=worker_script,
        pnfs_output_dir=pnfs_output_dir, local_output_dir=local_output_dir,
        inputs=inputs, outputs=outputs, extra=extra,
    )

    # Write command.json + the initial (pending) record *before* running, so a
    # slow/hanging submit still leaves a discoverable, jq-queryable artefact.
    command_file.write_text(json.dumps({
        "jobid": jobid, "runtype": runtype, "stem": stem,
        "command": cmd, "dry_run": dry_run, "n_jobs": n_jobs,
        "submitted": record["submitted"],
    }, indent=2) + "\n")
    records.atomic_write_json(gridlog, record)

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=timeout, env=build_submit_env(),
        )
        combined = (proc.stdout or "") + "\n--- stderr ---\n" + (proc.stderr or "")
        rc = proc.returncode
    except subprocess.TimeoutExpired as e:
        combined = (f"jobsub_submit timed out after {timeout}s\n"
                    f"stdout:\n{e.stdout or ''}\nstderr:\n{e.stderr or ''}")
        rc = -1
    except FileNotFoundError as e:
        combined = f"jobsub_submit not found: {e}\n"
        rc = -1
    submit_log.write_text(combined)

    cluster_id, status, finished, error = "", "submitted", None, ""
    if dry_run:
        status = "pending"
    elif rc != 0:
        status, finished = "failed", records.utc_now_iso()
        error = f"jobsub_submit exited {rc}"
    else:
        try:
            cluster_id = parse_cluster_id(combined)
        except ValueError as e:
            status, finished = "failed", records.utc_now_iso()
            error = str(e)

    updates = {"cluster_id": cluster_id, "status": status, "finished": finished}
    if error:
        updates["error"] = error
    record = records.update_record(gridlog, **updates)
    return record, gridlog
