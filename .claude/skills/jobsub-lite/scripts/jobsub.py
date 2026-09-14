#!/usr/bin/env python3
"""jobsub-lite skill CLI — FNAL grid submission with RCDS tarball publishing.

Run with **system python3** (never through pixi/conda). All state lives in the
project's `.jobsub/` directory (created by `init`); the skill directory itself
is read-only and copyable between projects.

    python3 jobsub.py init [--group dune]
    python3 jobsub.py tarball build --build-dir . --include sub1 sub2 [...]
    python3 jobsub.py tarball list
    python3 jobsub.py publish --tarball T.tar --label my_payload
    python3 jobsub.py verify --label my_payload
    python3 jobsub.py labels
    python3 jobsub.py adopt --label L --jobid J        # catalog an old publish
    python3 jobsub.py submit --worker w.sh -N 10 --tar-label my_payload \
        -- -R @TAR_DIR@ -O @PNFS_OUT@ <worker flags...>
    python3 jobsub.py status|cancel|fetchlog <jobid>
    python3 jobsub.py list [--active]
    python3 jobsub.py pull <jobid> [--suffix .root]

Placeholders substituted in worker args at submit time:
    @TAR_DIR@   the published payload's CVMFS dir (needs --tar-label)
    @PNFS_OUT@  the PNFS output dir (from --pnfs-out, or a default under
                <pnfs_scratch_base>/$USER/jobsub-lite/<runtype>/<stem>)
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

_SCRIPTS_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS_ROOT))

from lib import config, control, monitor, outputs, publish, records, tarball  # noqa: E402
from lib.submit import submit as do_submit                                    # noqa: E402

TAR_TOKEN = "@TAR_DIR@"
PNFS_TOKEN = "@PNFS_OUT@"


def _warn(msg: str) -> None:
    sys.stderr.write(f"warning: {msg}\n")


def _err(msg: str) -> int:
    sys.stderr.write(f"error: {msg}\n")
    return 2


# ── init ───────────────────────────────────────────────────────────────────────

def cmd_init(args) -> int:
    base = Path(args.dir or Path.cwd())
    state = config.ensure_state_dir(base)
    cfg_path = state / "config.json"

    if cfg_path.exists() and not args.force:
        print(f"state:  {state} (existing)")
        print(f"config: {cfg_path} (kept; --force to regenerate)")
    else:
        group = config.resolve_group(args.group)
        if not group:
            return _err("cannot determine the experiment group — pass --group "
                        "(e.g. --group dune) or set $JOBSUB_GROUP/$EXPERIMENT")
        cfg = config.default_config(group)
        if args.requirements:
            cfg["append_condor_requirements"] = args.requirements
        import json
        cfg_path.write_text(json.dumps(cfg, indent=2, sort_keys=True) + "\n")
        for key, val in cfg.items():
            if key.endswith("_bin") or key == "jobsub_bin":
                if not Path(val).is_absolute():
                    _warn(f"{key}: '{val}' not found on this host — edit {cfg_path}")
        print(f"state:  {state} (created)" if not args.force else f"state:  {state}")
        print(f"config: {cfg_path} (group={group})")

    # Auth report — informational only; jobsub_lite fetches a token on demand
    # but needs valid kerberos to do so.
    try:
        krb = subprocess.run(["klist", "-s"], capture_output=True).returncode == 0
    except FileNotFoundError:
        krb = False
    print(f"auth:   kerberos {'OK' if krb else 'MISSING — run kinit'}", end="")
    tok = os.environ.get("BEARER_TOKEN_FILE") or f"/run/user/{os.getuid()}/bt_u{os.getuid()}"
    print(f"; bearer token {'present' if Path(tok).exists() else 'absent (auto-fetched at submit)'}")
    return 0


# ── tarball ────────────────────────────────────────────────────────────────────

def cmd_tarball_build(args) -> int:
    res = tarball.build_tarball(
        build_dir=args.build_dir,
        toplevel_candidates=args.include,
        exclude_components=args.exclude_component or (),
        exclude_prefixes=args.exclude_prefix or (),
        exclude_suffixes=tuple(args.exclude_suffix or ()),
        name_prefix=args.name_prefix,
        output_path=args.output,
        force=args.force,
        background=args.background,
    )
    if "error" in res:
        return _err(res["error"])
    print(res["message"])
    print(f"tarball: {res['tarball_path']}")
    missing = set(args.include) - set(res["files_included"])
    if missing:
        _warn(f"not found under {args.build_dir} (skipped): {sorted(missing)}")
    return 0


def cmd_tarball_list(_args) -> int:
    d = tarball.tarball_dir()
    if not d.exists():
        print("(no tarballs)")
        return 0
    rows = sorted(d.glob("*.tar"), key=lambda p: p.stat().st_mtime, reverse=True)
    for p in rows:
        print(f"{p.stat().st_size / 1e6:10.1f} MB  {p.name}")
    print(f"\n{len(rows)} tarball(s) in {d}")
    return 0


# ── publish / catalog ──────────────────────────────────────────────────────────

def cmd_publish(args, cfg) -> int:
    tar_path = Path(args.tarball).resolve()
    if not tar_path.exists():
        return _err(f"tarball not found: {tar_path}")
    # Content hash makes same-label re-publishes decidable (no-op vs refuse)
    # *before* the sentinel job is submitted.
    local_sha = publish.file_sha256(tar_path)
    print(f"publishing {tar_path.name} to RCDS as '{args.label}' "
          f"(one sentinel grid job; typically 5–20 min)...")
    res = publish.publish_and_catalog(
        cfg, tarball_path=str(tar_path), label=args.label, local_sha=local_sha,
        size_mb=round(tar_path.stat().st_size / 1e6, 1),
        overwrite=args.overwrite, description=args.description or "",
        sentinel_script=args.sentinel,
    )
    if "error" in res:
        return _err(f"{res['error']}\n  publish log: {res.get('publish_log', '-')}")
    if res.get("message"):
        print(res["message"])
    print(f"label:     {res['label']}")
    print(f"cvmfs_dir: {res['cvmfs_dir']}")
    print(f"published: {res['published']}")
    return 0


def cmd_verify(args, _cfg) -> int:
    entry = publish.lookup_catalog(args.label)
    if entry is None:
        return _err(f"label not in catalog: {args.label!r} (see `labels`)")
    res = publish.verify_cvmfs(entry)
    cat = publish.load_catalog()
    cat["entries"][args.label] = entry
    publish.save_catalog(cat)
    print(f"label:          {args.label}")
    print(f"cvmfs_dir:      {entry.get('cvmfs_dir', '-')}")
    print(f"status:         {res['status']} ({res['reason']})")
    print(f"age:            {res['age_days']} day(s)")
    print(f"recommendation: {res['recommendation']}")
    return 0 if res["recommendation"] in ("ok", "warn") else 1


def cmd_labels(_args, _cfg) -> int:
    entries = publish.load_catalog()["entries"]
    if not entries:
        print("(catalog empty)")
        return 0
    for label, e in sorted(entries.items()):
        print(f"{label:<28} {e.get('published', '')[:10]:<12} {e.get('cvmfs_dir', '')}")
    print(f"\n{len(entries)} label(s)")
    return 0


def cmd_adopt(args, cfg) -> int:
    res = publish.label_from_job(cfg, label=args.label, jobid=args.jobid,
                                 description=args.description or "",
                                 overwrite=args.overwrite)
    if "error" in res:
        return _err(res["error"])
    print(f"label:     {res['label']}")
    print(f"cvmfs_dir: {res['cvmfs_dir']}")
    return 0


# ── submit ─────────────────────────────────────────────────────────────────────

def cmd_submit(args, cfg) -> int:
    worker = Path(args.worker).resolve()
    if not worker.exists():
        return _err(f"worker script not found: {worker}")
    if not os.access(worker, os.X_OK):
        return _err(f"worker script not executable: {worker} (chmod +x it)")

    now = datetime.now()
    stem = args.stem or f"{worker.stem}_{now.strftime('%Y%m%d-%H%M%S')}"
    user = os.environ.get("USER", "")

    worker_args = list(args.worker_args)
    if worker_args and worker_args[0] == "--":   # argparse keeps the separator
        worker_args = worker_args[1:]

    # Resolve the payload: a published catalog label, an explicit CVMFS dir, or
    # a local tar for dropbox.
    cvmfs_dir, tar_file_name = "", None
    if sum(bool(x) for x in (args.tar_label, args.tar_file, args.tar_dir)) > 1:
        return _err("--tar-label, --tar-dir and --tar-file are mutually exclusive")
    if args.tar_dir:
        cvmfs_dir = args.tar_dir.rstrip("/")
        if not Path(cvmfs_dir).is_dir():
            _warn(f"--tar-dir not visible on this host (may still exist on workers): {cvmfs_dir}")
        if not any(TAR_TOKEN in a for a in worker_args):
            _warn(f"--tar-dir given but no {TAR_TOKEN} placeholder in worker args — "
                  f"the worker will not receive the payload path")
    if args.tar_label:
        entry = publish.lookup_catalog(args.tar_label)
        if entry is None:
            return _err(f"label not in catalog: {args.tar_label!r} (see `labels`)")
        chk = publish.verify_cvmfs(entry)
        if chk["status"] == "missing":
            return _err(f"published payload for '{args.tar_label}' no longer on CVMFS "
                        f"(RCDS GC?) — re-publish before submitting")
        if chk["recommendation"] != "ok":
            _warn(f"payload '{args.tar_label}' age {chk['age_days']}d: {chk['recommendation']}")
        cvmfs_dir = entry["cvmfs_dir"]
        if not any(TAR_TOKEN in a for a in worker_args):
            _warn(f"--tar-label given but no {TAR_TOKEN} placeholder in worker args — "
                  f"the worker will not receive the payload path")
    if args.tar_file:
        tar_file_name = f"dropbox://{Path(args.tar_file).resolve()}"

    # PNFS output dir: substituted into worker args (and recorded) only when the
    # worker actually takes it, so records never point at never-written dirs.
    pnfs_out = ""
    uses_pnfs = args.pnfs_out or any(PNFS_TOKEN in a for a in worker_args)
    if uses_pnfs:
        base = cfg.get("pnfs_scratch_base", "")
        if args.pnfs_out:
            pnfs_out = args.pnfs_out
        elif base and user:
            pnfs_out = f"{base}/{user}/jobsub-lite/{args.runtype}/{stem}"
        else:
            return _err("cannot build a default PNFS output dir "
                        "(no pnfs_scratch_base in config / $USER unset) — pass --pnfs-out")
    if args.pnfs_out and not any(PNFS_TOKEN in a for a in worker_args):
        _warn(f"--pnfs-out given but no {PNFS_TOKEN} placeholder in worker args — "
              f"the worker will not receive the output dir")
    if any(TAR_TOKEN in a for a in worker_args) and not cvmfs_dir:
        _warn(f"{TAR_TOKEN} appears in worker args but no --tar-label/--tar-dir "
              f"was given — it will be substituted with an empty string")

    def _sub(a: str) -> str:
        return a.replace(TAR_TOKEN, cvmfs_dir).replace(PNFS_TOKEN, pnfs_out)
    worker_args = [_sub(a) for a in worker_args]

    group = args.group or cfg["default_group"]
    role = args.role or cfg.get("default_role", "Analysis")
    disk = args.disk or cfg.get("default_disk", "20GB")
    req = args.append_condor_requirements or cfg.get("append_condor_requirements")

    cmd = [cfg["jobsub_bin"], "-G", group, "--role", role,
           "--disk", disk, "-N", str(args.n_jobs)]
    if req:
        cmd += ["--append_condor_requirements", req]
    for f in (args.input_file or []):
        cmd += ["-f", f]
    if tar_file_name:
        cmd += ["--tar_file_name", tar_file_name]
    if args.memory:
        cmd += ["--memory", args.memory]
    if args.expected_lifetime:
        cmd += ["--expected-lifetime", args.expected_lifetime]
    for kv in (args.env or []):
        cmd += ["-e", kv]
    for extra_arg in (args.jobsub_arg or []):
        cmd.append(extra_arg)
    exe_index = len(cmd)
    cmd += [f"file://{worker}", *worker_args]

    extra = {"skill": "jobsub-lite"}
    if args.output_suffix:
        extra["output_suffix"] = args.output_suffix

    try:
        record, gridlog = do_submit(
            runtype=args.runtype, stem=stem, submit_cmd=cmd, n_jobs=args.n_jobs,
            submit_user=user, worker_script=str(worker),
            tarball_path=cvmfs_dir or (tar_file_name or ""),
            pnfs_output_dir=pnfs_out, extra=extra,
            dry_run=args.dry_run, exe_index=exe_index, when=now,
            timeout=int(cfg.get("submit_timeout_s", 600)),
        )
    except FileExistsError as e:
        return _err(str(e))

    print(f"jobid:  {record['jobid']}")
    line = f"status: {record['status']}"
    if record["cluster_id"]:
        line += f"  cluster: {record['cluster_id']}"
    print(line)
    print(f"record: {gridlog}")
    if pnfs_out:
        print(f"outputs: {pnfs_out}  (pull: jobsub.py pull {record['jobid']})")
    if record.get("error"):
        sys.stderr.write(f"error:  {record['error']}\n")
    return 1 if record["status"] == "failed" else 0


# ── job control ────────────────────────────────────────────────────────────────

def _fmt(rec: dict) -> str:
    return (f"{rec.get('jobid', ''):<46} {rec.get('status', ''):<10} "
            f"n={rec.get('n_jobs', 0):<5} done={rec.get('processes_done', 0):<5} "
            f"{rec.get('cluster_id', '') or '-'}")


def cmd_status(args, cfg) -> int:
    path = records.find_record_for_jobid(args.jobid)
    rec = monitor.refresh_status(path, cfg)
    print(f"jobid:      {rec['jobid']}")
    print(f"status:     {rec['status']}")
    print(f"cluster:    {rec.get('cluster_id') or '-'}")
    print(f"n_jobs:     {rec['n_jobs']}   processes_done: {rec.get('processes_done', 0)}"
          f" ({rec.get('processes_done_source') or '-'})")
    print(f"submitted:  {rec.get('submitted')}   finished: {rec.get('finished') or '-'}")
    print(f"record:     {path}")
    if rec.get("fetchlog_error"):
        print(f"fetchlog_error: {rec['fetchlog_error']}")
    cs = rec.get("cluster_state")
    if cs:
        print(f"queue:      total={cs.get('n_total', 0)} idle={cs.get('n_idle', 0)} "
              f"running={cs.get('n_running', 0)} held={cs.get('n_held', 0)} "
              f"done={cs.get('n_done', 0)}" + (f"  ({cs['error']})" if cs.get("error") else ""))
    return 0


def cmd_list(args, cfg) -> int:
    jobs = monitor.list_jobs(cfg, active_only=args.active)
    if not jobs:
        print("(no jobs)")
        return 0
    for rec in jobs:
        print(_fmt(rec))
    print(f"\n{len(jobs)} job(s)" + (" (active)" if args.active else ""))
    return 0


def cmd_cancel(args, cfg) -> int:
    path = records.find_record_for_jobid(args.jobid)
    rec = control.cancel(path, cfg)
    if not rec.get("_cancel_ok"):
        return _err(f"jobsub_rm failed for {rec['jobid']} — status kept as "
                    f"'{rec['status']}'. {rec.get('_cancel_detail', '')}")
    print(f"{rec['jobid']}: {rec['status']}  ({rec.get('_cancel_detail', '')})")
    return 0


def cmd_fetchlog(args, cfg) -> int:
    path = records.find_record_for_jobid(args.jobid)
    res = control.fetch_log(path, cfg)
    if "error" in res:
        return _err(res["error"])
    print(f"fetched {res['n_files']} file(s) -> {res['dest_dir']} (rc={res['returncode']})")
    return 0


def cmd_pull(args, cfg) -> int:
    path = records.find_record_for_jobid(args.jobid)
    # Resolve a drained-but-stale status first, so pull stamps (or declines to
    # stamp) against the campaign's real state.
    monitor.refresh_status(path, cfg)
    res = outputs.pull(path, cfg, suffix=args.suffix, overwrite=args.overwrite)
    if "error" in res:
        return _err(res["error"])
    print(res["message"])
    for f in res.get("failures", [])[:10]:
        sys.stderr.write(f"  FAIL: {f}\n")
    return 1 if res.get("status") == "failed" else 0


# ── main ───────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(
        description="jobsub-lite skill: FNAL grid submission + RCDS publishing")
    ap.add_argument("--config", default=None, help="path to config.json ($JOBSUB_CONFIG)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("init", help="create .jobsub/ + config in this project")
    p.add_argument("--dir", default=None, help="project root (default: CWD)")
    p.add_argument("--group", default=None, help="experiment group, e.g. dune")
    p.add_argument("--requirements", default=None,
                   help='append_condor_requirements, e.g. \'(TARGET.Microarch>="x86_64-v3")\'')
    p.add_argument("--force", action="store_true", help="regenerate config.json")

    p = sub.add_parser("tarball", help="build/list cached payload tarballs")
    tsub = p.add_subparsers(dest="tcmd", required=True)
    tb = tsub.add_parser("build", help="build (cache-keyed on tree mtimes)")
    tb.add_argument("--build-dir", required=True)
    tb.add_argument("--include", nargs="+", required=True,
                    help="top-level files/dirs (relative to --build-dir)")
    tb.add_argument("--exclude-component", action="append", metavar="NAME",
                    help="repeatable; skip any path containing this component (e.g. .git)")
    tb.add_argument("--exclude-prefix", action="append", metavar="PREFIX")
    tb.add_argument("--exclude-suffix", action="append", metavar="SUFFIX")
    tb.add_argument("--name-prefix", default="tarball")
    tb.add_argument("--output", default=None, help="explicit output .tar path")
    tb.add_argument("--force", action="store_true")
    tb.add_argument("--background", action="store_true",
                    help="detach the build (large payloads); watch the .log next to the .tar")
    tsub.add_parser("list", help="list cached tarballs")

    p = sub.add_parser("publish", help="upload a tarball to RCDS + catalog the CVMFS dir")
    p.add_argument("--tarball", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--description", default=None)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--sentinel", default=None, help="override sentinel worker script")

    p = sub.add_parser("verify", help="check a published label still exists / is fresh")
    p.add_argument("--label", required=True)

    sub.add_parser("labels", help="list the publish catalog")

    p = sub.add_parser("adopt", help="catalog an already-published tarball from a job's log")
    p.add_argument("--label", required=True)
    p.add_argument("--jobid", required=True)
    p.add_argument("--description", default=None)
    p.add_argument("--overwrite", action="store_true")

    p = sub.add_parser("submit", help="submit a worker script via jobsub_submit")
    p.add_argument("--worker", required=True)
    p.add_argument("-N", "--n-jobs", type=int, default=1)
    p.add_argument("--runtype", default="jobsub", help="record grouping label")
    p.add_argument("--stem", default=None, help="record stem (default: <worker>_<timestamp>)")
    p.add_argument("-G", "--group", default=None)
    p.add_argument("--role", default=None)
    p.add_argument("--disk", default=None, help="e.g. 20GB")
    p.add_argument("--memory", default=None, help="e.g. 2000MB")
    p.add_argument("--expected-lifetime", default=None, help="e.g. 8h")
    p.add_argument("--tar-label", default=None,
                   help=f"published payload label; its CVMFS dir replaces {TAR_TOKEN}")
    p.add_argument("--tar-file", default=None,
                   help="local .tar shipped via --tar_file_name dropbox:// (unpublished)")
    p.add_argument("--tar-dir", default=None,
                   help=f"explicit published CVMFS payload dir; replaces {TAR_TOKEN} "
                        f"(cross-project reuse without a catalog entry)")
    p.add_argument("-f", "--input-file", action="append", metavar="PATH",
                   help="repeatable; passed as `-f PATH` (use /pnfs paths or file://)")
    p.add_argument("--append-condor-requirements", default=None)
    p.add_argument("--env", action="append", metavar="KEY=VAL",
                   help="repeatable; passed as `-e KEY=VAL`")
    p.add_argument("--jobsub-arg", action="append", metavar="ARG",
                   help="repeatable; raw extra jobsub_submit option (escape hatch)")
    p.add_argument("--pnfs-out", default=None,
                   help=f"PNFS output dir substituted for {PNFS_TOKEN} in worker args")
    p.add_argument("--output-suffix", default=None,
                   help="expected output suffix (completion counting fallback)")
    p.add_argument("--dry-run", action="store_true", help="--no_submit; record as pending")
    p.add_argument("worker_args", nargs=argparse.REMAINDER,
                   help="args after `--` forwarded to the worker")

    p = sub.add_parser("status", help="re-poll + show one job"); p.add_argument("jobid")
    p = sub.add_parser("list", help="list jobs"); p.add_argument("--active", action="store_true")
    p = sub.add_parser("cancel", help="jobsub_rm a job"); p.add_argument("jobid")
    p = sub.add_parser("fetchlog", help="fetch worker logs"); p.add_argument("jobid")
    p = sub.add_parser("pull", help="ifdh-pull outputs to local disk")
    p.add_argument("jobid")
    p.add_argument("--suffix", default=None, help="only files ending with this")
    p.add_argument("--overwrite", action="store_true")

    args = ap.parse_args()

    if args.cmd == "init":
        return cmd_init(args)
    if args.cmd == "tarball":
        return cmd_tarball_build(args) if args.tcmd == "build" else cmd_tarball_list(args)

    cfg = config.load_config(args.config)
    return {
        "publish": cmd_publish, "verify": cmd_verify, "labels": cmd_labels,
        "adopt": cmd_adopt, "submit": cmd_submit, "status": cmd_status,
        "list": cmd_list, "cancel": cmd_cancel, "fetchlog": cmd_fetchlog,
        "pull": cmd_pull,
    }[args.cmd](args, cfg)


if __name__ == "__main__":
    sys.exit(main())
