"""Pull grid outputs from PNFS scratch back to local disk via ifdh.

Generic: *which* files to copy (`suffix`, defaulting to the record's
`extra.output_suffix`) and *how* to name them locally (`name_fn`) are
parameters. Walks `ifdh ls` of `pnfs_output_dir`'s per-process subdirs and
copies matching files into `local_output_dir`. Completion is counted in
*processes* (distinct subdirs with ≥1 matching file, mirroring
monitor.count_pnfs_outputs), and the record is stamped only for terminal,
failure-free pulls — see `pull`. All ifdh calls run under the scrubbed env.
(Vendored from genie-mcp, generalized.)
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Callable, Optional

from lib import records
from lib.submit_env import build_submit_env


def _ifdh_ls(path: str) -> tuple[bool, list[str]]:
    try:
        p = subprocess.run(["ifdh", "ls", path], capture_output=True,
                           text=True, timeout=120, env=build_submit_env())
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False, []
    if p.returncode != 0:
        return False, []
    return True, [l.strip() for l in (p.stdout or "").splitlines() if l.strip()]


def _ifdh_cp(src: str, dst: Path) -> tuple[bool, str]:
    try:
        p = subprocess.run(["ifdh", "cp", src, str(dst)], capture_output=True,
                           text=True, timeout=600, env=build_submit_env())
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        return False, str(e)
    return (p.returncode == 0), (p.stderr or "")[-200:]


def _default_name(process_str: str, basename: str) -> str:
    """Prefix with the process subdir so files from different processes don't
    collide in the flat local dir."""
    return f"{process_str}_{basename}" if process_str else basename


def pull(
    record_path: Path,
    cfg: dict,
    *,
    suffix: Optional[str] = None,
    name_fn: Optional[Callable[[str, str], str]] = None,
    overwrite: bool = False,
) -> dict:
    """Copy outputs from `pnfs_output_dir` into `local_output_dir`.

    `suffix`  — only copy files ending with it (None = the record's
                `extra.output_suffix` if set, else every file).
    `name_fn` — (process_str, basename) -> local filename (default: prefix the
                process subdir name).

    Stamping contract: `outputs_pulled` / `processes_done` / `status` are
    written only when the record is already terminal **and** no copy failed —
    a pull on an in-flight campaign copies what exists and records
    `local_output_dir`, but must not freeze the monitoring verdict
    (refresh_status short-circuits terminal states), and a failed transfer
    must stay resumable without `--overwrite`.
    """
    record = records.read_record(record_path)
    name_fn = name_fn or _default_name
    if suffix is None:
        suffix = record.get("extra", {}).get("output_suffix") or None
    terminal = record.get("status") in records.TERMINAL_STATES

    if record.get("outputs_pulled") and not overwrite:
        return {
            "status": "already_pulled", "job_id": record["jobid"],
            "local_output_dir": record.get("local_output_dir", ""),
            "n_pulled": record.get("processes_done", 0), "n_expected": record["n_jobs"],
            "message": "outputs already pulled; pass overwrite=True to re-copy",
        }

    pnfs_dir = record.get("pnfs_output_dir", "")
    if not pnfs_dir:
        return {"error": "record has no pnfs_output_dir", "job_id": record["jobid"]}

    local_dir = Path(record.get("local_output_dir")
                     or (Path(record_path).parent / f"{record['stem']}.outputs"))
    local_dir.mkdir(parents=True, exist_ok=True)

    ok, top = _ifdh_ls(pnfs_dir)
    if not ok:
        return {"error": f"ifdh ls failed for {pnfs_dir} (does it exist? is ifdh on PATH?)",
                "job_id": record["jobid"]}

    process_dirs = [p for p in top if p.endswith("/") and p.rstrip("/") != pnfs_dir.rstrip("/")]
    pulled: list[str] = []
    failures: list[str] = []
    procs_with_files: set[str] = set()

    for pdir in process_dirs:
        process_str = pdir.rstrip("/").rsplit("/", 1)[-1]
        ok2, files = _ifdh_ls(pdir)
        if not ok2:
            failures.append(f"ifdh ls failed: {pdir}")
            continue
        for fpath in files:
            if fpath.endswith("/"):
                continue
            basename = fpath.rsplit("/", 1)[-1]
            if suffix and not basename.endswith(suffix):
                continue
            local_path = local_dir / name_fn(process_str, basename)
            if local_path.exists() and not overwrite:
                pulled.append(str(local_path))
                procs_with_files.add(process_str)
                continue
            cok, detail = _ifdh_cp(fpath, local_path)
            if cok:
                pulled.append(str(local_path))
                procs_with_files.add(process_str)
            else:
                failures.append(f"ifdh cp {fpath}: {detail}")

    n_files = len(pulled)
    n_procs = len(procs_with_files)
    n_expected = record["n_jobs"]

    updates: dict = {"local_output_dir": str(local_dir)}
    if terminal and not failures:
        updates["outputs_pulled"] = True
        updates["processes_done"] = n_procs
        if record.get("status") != "cancelled":
            updates["status"] = ("done" if n_procs == n_expected
                                 else ("partial" if n_procs > 0 else "failed"))
    records.update_record(record_path, **updates)

    status = updates.get("status", record.get("status"))
    if not terminal:
        msg = (f"campaign still active (status={status}) — copied {n_files} file(s) "
               f"from {n_procs} process(es); record not stamped, re-run pull after "
               f"the queue drains")
    else:
        msg = (f"pulled {n_files} file(s) across {n_procs}/{n_expected} process(es) "
               f"into {local_dir}")
        if failures:
            msg += (f"; {len(failures)} failure(s) — record not stamped, "
                    f"re-run pull to resume")
    return {
        "job_id": record["jobid"], "n_pulled": n_files, "n_processes": n_procs,
        "n_expected": n_expected, "local_output_dir": str(local_dir),
        "files": pulled, "failures": failures, "status": status, "message": msg,
    }
