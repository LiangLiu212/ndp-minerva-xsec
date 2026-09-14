"""Publish a tarball to RCDS/CVMFS and track it in a label->entry catalog.

RCDS assigns each upload to a random /cvmfs/fifeuserN.opensciencegrid.org/sw/
<group>/<hash>/ (N=1..4), so we recover the real path by running one **sentinel
grid job** (templates/publish_only.sh) that echoes
`PUBLISH_SENTINEL_CVMFS_DIR=<path>`, then fetch its log and grep that out. The
catalog (`<state_dir>/catalog.json`) maps a label to the published CVMFS path +
a publish timestamp; `verify_cvmfs` flags staleness (RCDS garbage-collects
~30d). The experiment group is read from the config, not hard-coded.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from lib import config, monitor, records
from lib.submit import parse_cluster_id
from lib.submit_env import build_submit_env

_SKILL_ROOT = Path(__file__).resolve().parents[2]
_DEFAULT_SENTINEL = _SKILL_ROOT / "templates" / "publish_only.sh"

# Require an absolute path so we match the sentinel's *stdout* (the expanded
# /cvmfs/... path), not the fetched copy of the worker script's own source line
# (`PUBLISH_SENTINEL_CVMFS_DIR=${INPUT_TAR_FILE%/}`).
_RE_SENTINEL_CVMFS = re.compile(r"PUBLISH_SENTINEL_CVMFS_DIR=(/\S+)")
_CVMFS_RCDS_REPOS = tuple(
    f"/cvmfs/fifeuser{i}.opensciencegrid.org" for i in (1, 2, 3, 4)
)
_GC_WARN_DAYS = 21
_GC_FAIL_DAYS = 28


def catalog_path() -> Path:
    return config.state_dir() / "catalog.json"


def publish_logs_dir() -> Path:
    return config.state_dir() / "runs" / "_publish"


def cvmfs_roots(group: str) -> tuple:
    return tuple(f"{repo}/sw/{group}" for repo in _CVMFS_RCDS_REPOS)


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_rcds_hash(text: str, group: str) -> Optional[str]:
    m = re.findall(rf"Publishing hash {re.escape(group)}/([0-9a-f]+)", text)
    return m[-1] if m else None


def file_sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ── catalog (label -> entry) ──────────────────────────────────────────────────

def load_catalog(path: Path | None = None) -> dict:
    path = path or catalog_path()
    if not Path(path).exists():
        return {"entries": {}}
    try:
        data = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError):
        return {"entries": {}}
    if not isinstance(data.get("entries"), dict):
        data["entries"] = {}
    return data


def save_catalog(catalog: dict, path: Path | None = None) -> None:
    path = Path(path or catalog_path())
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".catalog.", suffix=".json", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(catalog, f, indent=2, sort_keys=True)
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def lookup_catalog(label: str, path: Path | None = None) -> Optional[dict]:
    return load_catalog(path)["entries"].get(label)


def add_to_catalog(entry: dict, *, overwrite: bool = False, path: Path | None = None) -> dict:
    catalog = load_catalog(path)
    label = entry["label"]
    existing = catalog["entries"].get(label)
    if existing is not None and not overwrite:
        # Identical *known* content is an idempotent no-op. An empty sha means
        # "content unknown" and must never match — otherwise every same-label
        # re-publish would silently keep the old entry.
        if entry.get("local_sha") and existing.get("local_sha") == entry.get("local_sha"):
            return existing
        raise ValueError(
            f"label '{label}' already in catalog with sha {existing.get('local_sha')!r}; "
            f"pass overwrite=True to replace"
        )
    catalog["entries"][label] = entry
    save_catalog(catalog, path)
    return entry


# ── verify / staleness ──────────────────────────────────────────────────────

def _parse_iso(s: str) -> Optional[datetime]:
    if not s:
        return None
    try:
        dt = datetime.fromisoformat(s)
    except ValueError:
        return None
    return dt.replace(tzinfo=timezone.utc) if dt.tzinfo is None else dt


def verify_cvmfs(entry: dict) -> dict:
    """Best-effort existence + age check for a catalog entry's cvmfs_dir."""
    cvmfs_dir = entry.get("cvmfs_dir", "")
    published = _parse_iso(entry.get("published", ""))
    now = datetime.now(timezone.utc)
    age_days = (now - published).days if published else -1

    # RCDS unpacks the tarball on upload, so the published <hash>/ dir holds the
    # extracted tree — there is no .tar on CVMFS. Verify the dir exists.
    if cvmfs_dir and Path(cvmfs_dir).is_dir():
        status, reason = "exists", "stat ok"
    elif not any(Path(r).exists() for r in _CVMFS_RCDS_REPOS):
        status, reason = "unknown", "no CVMFS fifeuserN repo mounted on this host"
    else:
        status, reason = "missing", f"path not present in CVMFS: {cvmfs_dir}"

    if age_days < 0:
        rec = "warn"
    elif age_days > _GC_FAIL_DAYS or status == "missing":
        rec = "republish"
    elif age_days > _GC_WARN_DAYS:
        rec = "warn"
    else:
        rec = "ok"

    result = {"status": status, "reason": reason, "age_days": age_days, "recommendation": rec}
    entry["last_verified"] = _now_iso()
    entry["last_verified_result"] = result
    return result


# ── publish via sentinel grid job ──────────────────────────────────────────────

def grep_sentinel_cvmfs(fetch_dir: Path) -> Optional[str]:
    """Return the PUBLISH_SENTINEL_CVMFS_DIR path from any fetched log, or None."""
    for f in sorted(Path(fetch_dir).rglob("*")):
        if not f.is_file():
            continue
        m = _RE_SENTINEL_CVMFS.search(f.read_text(errors="replace"))
        if m:
            return m.group(1).rstrip("/")
    return None


def locate_cvmfs_by_hash(rcds_hash: str, group: str) -> Optional[str]:
    """RCDS unpacks each upload to /cvmfs/fifeuserN/sw/<group>/<hash>/. Return
    that dir if present on this host, else None — a fetchlog-independent
    fallback."""
    if not rcds_hash:
        return None
    for root in cvmfs_roots(group):
        cand = f"{root}/{rcds_hash}"
        if Path(cand).is_dir():
            return cand
    return None


def publish_to_cvmfs(
    cfg: dict,
    *,
    tarball_path: str,
    label: str,
    sentinel_script: Optional[str] = None,
    log_prefix: str = "",
    poll_timeout_s: int = 1800,
    poll_interval_s: int = 30,
    fetch_attempts: int = 3,
    fetch_interval_s: int = 30,
) -> dict:
    """Upload `tarball_path` to RCDS by running one sentinel grid job; return
    {rcds_hash, cluster_id, cvmfs_dir, publish_log, fetched_dir, command_str}
    or {error, ...}."""
    group = cfg["default_group"]
    sentinel = Path(sentinel_script or _DEFAULT_SENTINEL)
    if not sentinel.exists():
        return {"error": f"sentinel worker not found: {sentinel}"}
    if not os.access(sentinel, os.X_OK):
        return {"error": f"sentinel worker not executable: {sentinel}"}

    # Short lifetime class: the sentinel only echoes env vars, and shorter
    # classes schedule faster.
    cmd = [cfg["jobsub_bin"], "-G", group,
           "--role", cfg.get("default_role", "Analysis"), "-N", "1",
           "--expected-lifetime", "1h",
           "--tar_file_name", f"dropbox://{tarball_path}", f"file://{sentinel}"]

    logs_dir = publish_logs_dir()
    logs_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = logs_dir / f"{log_prefix}{label}_{ts}.log"

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=900, env=build_submit_env())
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        log_path.write_text(f"jobsub_submit failed: {e}\n")
        return {"error": f"jobsub_submit failed: {e}", "publish_log": str(log_path)}

    combined = (proc.stdout or "") + "\n--- stderr ---\n" + (proc.stderr or "")
    log_path.write_text(combined)
    if proc.returncode != 0:
        return {"error": f"jobsub_submit exited {proc.returncode}",
                "publish_log": str(log_path), "command_str": " ".join(cmd)}

    rcds_hash = parse_rcds_hash(combined, group)
    if not rcds_hash:
        return {"error": "RCDS publish hash not found in jobsub_submit output",
                "publish_log": str(log_path), "command_str": " ".join(cmd)}
    try:
        cluster_id = parse_cluster_id(combined)
    except ValueError:
        return {"error": "cluster id not found in jobsub_submit output",
                "publish_log": str(log_path), "rcds_hash": rcds_hash}

    # Poll until the sentinel cluster leaves the queue. A just-submitted job is
    # not yet registered with the schedd, so jobsub_q reports `empty` for a few
    # seconds — treat `empty` as "drained" only after we've seen the job appear
    # (or after an appearance grace window, in case it ran and drained between
    # polls).
    deadline = time.monotonic() + poll_timeout_s
    appear_deadline = time.monotonic() + max(poll_interval_s * 4, 120)
    seen = False
    while True:
        q = monitor.query_jobsub_status(cluster_id, cfg)
        # Trust only healthy polls (same guard as monitor.refresh_status): a
        # crashed or timed-out jobsub_q yields an empty aggregate that is
        # indistinguishable from a drained queue.
        healthy = not q.get("error") and q.get("raw_returncode", 0) == 0
        if healthy and q.get("empty"):
            if seen or time.monotonic() > appear_deadline:
                break
        elif healthy:
            seen = True
        if time.monotonic() > deadline:
            return {"error": f"sentinel {cluster_id} did not leave queue in {poll_timeout_s}s",
                    "publish_log": str(log_path), "rcds_hash": rcds_hash, "cluster_id": cluster_id}
        time.sleep(poll_interval_s)

    # Fetch worker log + grep the real CVMFS dir. The sentinel's logs are not
    # always retrievable the instant the job drains, so retry the fetch+grep a
    # few times; and since we already parsed the RCDS hash, fall back to locating
    # the unpacked /cvmfs/.../<hash> dir directly (fetchlog-independent).
    fetch_dir = logs_dir / f"{log_prefix}{label}_{ts}_fetched"
    fetch_dir.mkdir(parents=True, exist_ok=True)
    fetch_cmd = [cfg["jobsub_fetchlog_bin"], "--jobid", cluster_id,
                 "-G", group, "--unzipdir", str(fetch_dir)]

    cvmfs_dir = None
    fetch_err = ""
    for attempt in range(max(1, fetch_attempts)):
        if attempt:
            time.sleep(fetch_interval_s)
        try:
            subprocess.run(fetch_cmd, capture_output=True, text=True,
                           timeout=600, env=build_submit_env())
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            fetch_err = f"jobsub_fetchlog failed: {e}"
            continue
        cvmfs_dir = grep_sentinel_cvmfs(fetch_dir)
        if cvmfs_dir:
            break

    if not cvmfs_dir:
        # Fallback: the RCDS publish already happened; find the unpacked dir.
        cvmfs_dir = locate_cvmfs_by_hash(rcds_hash, group)

    if not cvmfs_dir:
        detail = f" ({fetch_err})" if fetch_err else ""
        return {"error": f"PUBLISH_SENTINEL_CVMFS_DIR not found in fetched logs under "
                         f"{fetch_dir} after {max(1, fetch_attempts)} attempt(s){detail}, and "
                         f"/cvmfs/fifeuser{{1..4}}/sw/{group}/{rcds_hash} not visible on this host yet",
                "publish_log": str(log_path), "rcds_hash": rcds_hash,
                "cluster_id": cluster_id, "fetched_dir": str(fetch_dir)}

    return {"rcds_hash": rcds_hash, "cluster_id": cluster_id, "cvmfs_dir": cvmfs_dir,
            "group": group, "publish_log": str(log_path), "fetched_dir": str(fetch_dir),
            "command_str": " ".join(cmd)}


def publish_and_catalog(
    cfg: dict,
    *,
    tarball_path: str,
    label: str,
    local_sha: str = "",
    size_mb: Optional[float] = None,
    overwrite: bool = False,
    description: str = "",
    sentinel_script: Optional[str] = None,
    log_prefix: str = "",
    extra: Optional[dict] = None,
) -> dict:
    """Publish `tarball_path` then record the CVMFS path in the catalog under
    `label`. Returns the catalog entry, or {error}."""
    if not label or "/" in label or label.startswith("."):
        return {"error": f"invalid label: {label!r}"}

    # Guard the label *before* the expensive sentinel job: a same-label publish
    # without overwrite must either be a provable no-op (identical content,
    # still on CVMFS) or refuse now — not after 5–20 min of grid time.
    existing = lookup_catalog(label)
    refresh_same_content = False
    if existing is not None and not overwrite:
        same = bool(local_sha) and existing.get("local_sha") == local_sha
        if not same:
            return {"error": (
                f"label '{label}' already in catalog "
                f"(published {existing.get('published', '?')[:10]}, "
                f"sha {existing.get('local_sha') or 'unrecorded'}); the new tarball's "
                f"content differs or is unhashed — pass --overwrite to replace")}
        chk = verify_cvmfs(existing)
        if chk["status"] == "exists":
            return {**existing, "no_op": True,
                    "message": "identical content already published — sentinel skipped (no-op)"}
        if chk["status"] == "unknown":
            return {**existing, "no_op": True,
                    "message": ("identical content already in the catalog; CVMFS is not "
                                "mounted here so the published dir cannot be verified — "
                                "pass --overwrite to force a re-publish")}
        refresh_same_content = True  # same content but the dir is gone (RCDS GC): re-publish

    pub = publish_to_cvmfs(cfg, tarball_path=tarball_path, label=label,
                           sentinel_script=sentinel_script, log_prefix=log_prefix)
    if "error" in pub:
        return pub

    entry = {
        "label": label, "local_path": tarball_path, "local_sha": local_sha,
        "size_mb": size_mb, "rcds_hash": pub["rcds_hash"], "cluster_id": pub.get("cluster_id", ""),
        "cvmfs_dir": pub["cvmfs_dir"], "group": pub.get("group", cfg["default_group"]),
        "published": _now_iso(), "last_verified": "", "last_verified_result": {},
        "description": description, "publish_log": pub["publish_log"],
        "fetched_dir": pub.get("fetched_dir", ""), **(extra or {}),
    }
    try:
        return add_to_catalog(entry, overwrite=overwrite or refresh_same_content)
    except ValueError as e:
        return {"error": str(e), "publish_log": pub["publish_log"]}


def label_from_job(cfg: dict, *, label: str, jobid: str,
                   description: str = "", overwrite: bool = False) -> dict:
    """Adopt an already-published tarball into the catalog by parsing the RCDS
    hash from an existing job's submit log (needs the CVMFS path visible here)."""
    if not label or "/" in label or label.startswith("."):
        return {"error": f"invalid label: {label!r}"}
    group = cfg["default_group"]
    try:
        rec_path = records.find_record_for_jobid(jobid)
    except (FileNotFoundError, ValueError) as e:
        return {"error": str(e)}
    record = records.read_record(rec_path)
    submit_log = record.get("submit_log_file", "")
    if not submit_log or not Path(submit_log).exists():
        return {"error": f"submit log not found for job '{jobid}': {submit_log!r}"}
    rcds_hash = parse_rcds_hash(Path(submit_log).read_text(), group)
    if not rcds_hash:
        return {"error": f"no 'Publishing hash {group}/<hex>' line in submit log for '{jobid}'"}
    tarball_path = record.get("tarball_path", "")
    # RCDS unpacks the upload, so CVMFS holds the extracted <hash>/ dir, not the
    # .tar — locate the dir, not a file inside it.
    cvmfs_dir = locate_cvmfs_by_hash(rcds_hash, group)
    if not cvmfs_dir:
        return {"error": f"could not locate fifeuser{{1..4}}/sw/{group}/{rcds_hash} "
                         f"on this host; re-publish via publish_and_catalog instead"}
    entry = {
        "label": label, "local_path": tarball_path, "local_sha": "",
        "rcds_hash": rcds_hash, "cvmfs_dir": cvmfs_dir, "group": group,
        "published": record.get("submitted", _now_iso()), "last_verified": "",
        "last_verified_result": {}, "description": description, "adopted_from_job": jobid,
    }
    try:
        return add_to_catalog(entry, overwrite=overwrite)
    except ValueError as e:
        return {"error": str(e)}
