"""Campaign bookkeeping for grid processing of the MINERvA Open Data AnaTuples.

A campaign is a directory `grid/campaigns/<name>/` with

    campaign.json                    every file: tag, url, beam, kind, playlist, harvested POT, status,
                                     worklist, cluster/process, attempts
    worklists/<beam>_<kind>_<pl>.txt one URL per line, the `-f` input of a job cluster
                                     (`worker.sh -W` slices it by $PROCESS and files-per-process)
    resubmit_<n>.txt                 worklists of files that are missing / failed / incomplete

Commands (python -m ndp grid ...): harvest-pot, plan, status, resubmit, harvest, merge.
Status/harvest read the workers' outputs on PNFS through `xrdfs`/`xrdcp` against the DUNE
dCache door (a bearer token from `htgettoken` must be present) or, when an `ifdh` is on
PATH, through ifdh. Workers never consult PNFS: the worklist is the source of truth.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from ..config import REPO_ROOT

FILELISTS = REPO_ROOT / ".claude/skills/minerva-open-data/scripts/filelists"
POT_TABLE = REPO_ROOT / "resources/minerva/opendata_pot_per_file.tsv"
CAMPAIGNS = REPO_ROOT / "grid/campaigns"
DUNE_DOOR = "root://fndca1.fnal.gov:1094"
#: the system XRootD client carries the SciTokens auth plugin; the conda one in a pixi env may not
XRDFS = "/usr/bin/xrdfs" if os.path.exists("/usr/bin/xrdfs") else "xrdfs"
XRDCP = "/usr/bin/xrdcp" if os.path.exists("/usr/bin/xrdcp") else "xrdcp"
_LIST_RE = re.compile(r"MediumEnergy_(FHC|RHC)_(Data|StandardMC)_Playlist(\w+?)\.txt$")


# ---- file lists and POT --------------------------------------------------------------------------
def list_files(beams=("FHC",), kinds=("Data", "StandardMC"), playlists=None) -> list[dict]:
    """Every (beam, kind, playlist, url) of the published file lists (downloaded by minerva_od.py filelists)."""
    if not FILELISTS.exists():
        raise FileNotFoundError(f"{FILELISTS} missing: run `python3 .claude/skills/minerva-open-data/scripts/minerva_od.py filelists`")
    from ..adapters.minerva_anatuple import cache_tag
    out = []
    for f in sorted(FILELISTS.glob("*.txt")):
        m = _LIST_RE.search(f.name)
        if not m:
            continue
        beam, kind, pl = m.groups()
        if pl.startswith("6I_") or beam not in beams or kind not in kinds or (playlists and pl not in playlists):
            continue
        for line in f.read_text().splitlines():
            line = line.strip()
            if line.startswith("root://"):
                out.append({"beam": beam, "kind": "mc" if kind == "StandardMC" else "data", "playlist": pl,
                            "url": line, "tag": cache_tag(line)})
    return out


def harvest_pot(out: Path = POT_TABLE, workers: int = 24, log=print) -> Path:
    """POT_Used / POT_Total / entry totals of every published data + StandardMC file (Meta tree, streamed)."""
    from ..adapters.minerva_anatuple import read_pot
    files = list_files(beams=("FHC", "RHC"))
    log(f"reading the Meta tree of {len(files)} files with {workers} threads")
    rows = []

    def one(item, tries=3):
        err = None
        for _ in range(tries):
            try:
                return {**item, **read_pot(item["url"]), "error": None}
            except Exception as e:  # noqa: BLE001
                err = f"{type(e).__name__}: {str(e)[:80]}"; time.sleep(2)
        return {**item, "pot_used": 0.0, "pot_total": 0.0, "n_meta_entries": 0, "error": err}

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for i, fu in enumerate(as_completed([ex.submit(one, it) for it in files])):
            rows.append(fu.result())
            if (i + 1) % 500 == 0:
                log(f"  {i + 1}/{len(files)} in {time.time() - t0:.0f} s")
    rows.sort(key=lambda r: (r["beam"], r["kind"], r["playlist"], r["tag"]))
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["beam", "kind", "playlist", "tag", "url", "pot_used", "pot_total", "n_meta_entries", "total_reco_entries", "total_truth_entries", "error"]
    with open(out, "w") as fh:
        fh.write("# MINERvA Open Data per-file POT, read from each file's Meta tree over XRootD on " + time.strftime("%Y-%m-%d") + "\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    n_fail = sum(1 for r in rows if r["error"])
    log(f"wrote {out}: {len(rows)} files, {n_fail} failures, {time.time() - t0:.0f} s")
    return out


def load_pot_table(path: Path = POT_TABLE) -> dict[str, dict]:
    if not path.exists():
        return {}
    rows = {}
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts; continue
            r = dict(zip(header, parts))
            rows[r["tag"]] = r
    return rows


# ---- campaign.json -----------------------------------------------------------------------------
def campaign_dir(name: str) -> Path:
    return CAMPAIGNS / name


def load_campaign(name: str) -> dict:
    p = campaign_dir(name) / "campaign.json"
    if not p.exists():
        raise FileNotFoundError(f"no campaign {name!r} ({p}); run `ndp grid plan` first")
    return json.loads(p.read_text())


def save_campaign(c: dict) -> Path:
    d = campaign_dir(c["name"]); d.mkdir(parents=True, exist_ok=True)
    p = d / "campaign.json"
    p.write_text(json.dumps(c, indent=1))
    return p


def plan(name: str, channel: str, beams=("FHC",), kinds=("data", "mc"), playlists=None, files_per_process=None,
         pnfs_base: str | None = None, log=print) -> dict:
    """Create the campaign: file table (with harvested POT) + one worklist per (beam, kind, playlist)."""
    from ..io import timestamp
    fpp = {"mc": 8, "data": 100}      # measured 2026-09-14: ~70 s per streamed MC file, ~4 s per data file -> ~10 / ~7 min per process
    fpp.update(files_per_process or {})
    kinds_up = tuple({"data": "Data", "mc": "StandardMC"}[k] for k in kinds)
    files = list_files(beams=beams, kinds=kinds_up, playlists=playlists)
    pot = load_pot_table()
    user = os.environ.get("USER", "xxx")
    c = {"name": name, "created": timestamp(), "channel": channel, "beams": list(beams), "kinds": list(kinds),
         "playlists": sorted({f["playlist"] for f in files}), "files_per_process": fpp,
         "pnfs_base": pnfs_base or f"/pnfs/dune/scratch/users/{user}/ndp-stream/{name}", "worklists": {}, "files": {}}
    d = campaign_dir(name); (d / "worklists").mkdir(parents=True, exist_ok=True)
    groups: dict[tuple, list] = {}
    for f in files:
        groups.setdefault((f["beam"], f["kind"], f["playlist"]), []).append(f)
    for (beam, kind, pl), items in sorted(groups.items()):
        wl = f"{beam}_{kind}_{pl}"
        p = d / "worklists" / f"{wl}.txt"
        p.write_text("".join(f["url"] + "\n" for f in items))
        n_proc = -(-len(items) // fpp[kind])
        c["worklists"][wl] = {"file": str(p), "beam": beam, "kind": kind, "playlist": pl, "n_files": len(items),
                              "files_per_process": fpp[kind], "n_processes": n_proc,
                              "pnfs_out": f"{c['pnfs_base']}/{wl}", "jobs": []}
        for i, f in enumerate(items):
            c["files"][f["tag"]] = {**f, "worklist": wl, "line": i, "process": i // fpp[kind], "status": "planned",
                                    "pot_used_harvest": float(pot.get(f["tag"], {}).get("pot_used") or 0.0) or None,
                                    "attempts": 0}
    save_campaign(c)
    log(f"campaign {name}: {len(files)} files in {len(groups)} worklists -> {d}")
    for wl, w in c["worklists"].items():
        log(f"  {wl}: {w['n_files']} files, {w['n_processes']} processes x {w['files_per_process']} files")
    return c


# ---- PNFS access (login node) ---------------------------------------------------------------------
def _have(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def xrootd_path(pnfs_path: str) -> str:
    """dCache's XRootD namespace: /pnfs/<exp>/... is served as /pnfs/fs/usr/<exp>/... on the door."""
    p = pnfs_path
    if p.startswith("/pnfs/fs/usr/"):
        return p
    if p.startswith("/pnfs/"):
        return "/pnfs/fs/usr/" + p[len("/pnfs/"):]
    return p


def _token_env() -> dict:
    """Environment for xrdfs/xrdcp: point BEARER_TOKEN_FILE at jobsub_lite's token when nothing else is set."""
    env = dict(os.environ)
    if not env.get("BEARER_TOKEN_FILE") or not Path(env["BEARER_TOKEN_FILE"]).exists():   # a stale variable hides a valid token
        env.pop("BEARER_TOKEN_FILE", None)
        for cand in (f"/tmp/bt_token_dune_Analysis_{os.getuid()}", f"/run/user/{os.getuid()}/bt_u{os.getuid()}", f"/tmp/bt_u{os.getuid()}"):
            if Path(cand).exists():
                env["BEARER_TOKEN_FILE"] = cand; break
    return env


def pnfs_ls(path: str) -> list[str]:
    """Names under a /pnfs path: ifdh if available, else xrdfs against the DUNE door (bearer token)."""
    if _have("ifdh"):
        r = subprocess.run(["ifdh", "ls", path, "1"], capture_output=True, text=True, timeout=300, env=_token_env())
        if r.returncode != 0:
            raise RuntimeError(f"ifdh ls {path}: {r.stderr.strip()[:200]}")
        return [Path(l.strip()).name for l in r.stdout.splitlines() if l.strip() and not l.strip().endswith(path.rstrip("/"))]
    r = subprocess.run([XRDFS, DUNE_DOOR, "ls", xrootd_path(path)], capture_output=True, text=True, timeout=300, env=_token_env())
    if r.returncode != 0:
        raise RuntimeError(f"xrdfs ls {path}: {r.stderr.strip()[:200]} (token: htgettoken -a htvaultprod.fnal.gov -i dune)")
    return [Path(l.strip()).name for l in r.stdout.splitlines() if l.strip()]


def pnfs_copy(src_pnfs: str, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if _have("ifdh"):
        r = subprocess.run(["ifdh", "cp", src_pnfs, str(dst)], capture_output=True, text=True, timeout=1800, env=_token_env())
    else:
        r = subprocess.run([XRDCP, "-f", "-s", f"{DUNE_DOOR}/{xrootd_path(src_pnfs)}", str(dst)], capture_output=True, text=True, timeout=1800, env=_token_env())
    if r.returncode != 0:
        raise RuntimeError(f"copy {src_pnfs}: {r.stderr.strip()[:200]}")


def scan_outputs(pnfs_out: str) -> dict[str, dict]:
    """{tag: {"dir": <pnfs process dir>, "files": [...]}} for every manifest_<tag>.json under a worklist's output base."""
    found = {}
    for proc in pnfs_ls(pnfs_out):
        if not re.fullmatch(r"\d{4}", proc):
            continue
        names = pnfs_ls(f"{pnfs_out}/{proc}")
        for n in names:
            m = re.fullmatch(r"manifest_(\w+)\.json", n)
            if m:
                found[m.group(1)] = {"dir": f"{pnfs_out}/{proc}", "files": names}
    return found


def status(name: str, pot_check: bool = True, log=print) -> dict:
    """Refresh every file's status from the sidecars on PNFS; returns counts per status."""
    c = load_campaign(name)
    d = campaign_dir(name); side_dir = d / "sidecars"; side_dir.mkdir(exist_ok=True)
    for wl, w in c["worklists"].items():
        mine = [f for f in c["files"].values() if f["worklist"] == wl]
        if all(f["status"] == "planned" for f in mine) or all(f["status"] == "done" for f in mine):
            continue                                       # nothing submitted yet, or already complete
        try:
            found = scan_outputs(w["pnfs_out"])
        except RuntimeError as e:
            log(f"{wl}: cannot list {w['pnfs_out']}: {e}")
            continue
        for tag, info in found.items():
            if tag not in c["files"]:
                continue
            f = c["files"][tag]
            local = side_dir / f"manifest_{tag}.json"
            if not local.exists():
                pnfs_copy(f"{info['dir']}/manifest_{tag}.json", local)
            s = json.loads(local.read_text())
            f["pnfs_dir"] = info["dir"]; f["attempts"] = int(s.get("attempts", f.get("attempts", 0)) or 0)
            f["cluster"], f["process_id"] = s.get("cluster"), s.get("process")
            if s.get("status") != "ok":
                f["status"] = "failed"; f["error"] = (s.get("errors") or [""])[-1]
                continue
            problems = []
            expect = [f"reco_{tag}.npz"] + ([f"truth_{tag}.npz", f"reco_{tag}_truthcols.npz", f"truth_{tag}_skim.npz", f"reco_{tag}_truthcols_skim.npz"] if f["kind"] == "mc" else [])
            missing = [n for n in expect if n not in info["files"]]
            if missing:
                problems.append(f"missing outputs {missing}")
            if s.get("total_reco_entries") is not None and s.get("n_reco") is not None and int(s["total_reco_entries"]) != int(s["n_reco"]):
                problems.append(f"n_reco {s['n_reco']} != Meta Total_Reco_Entries {s['total_reco_entries']}")
            if s.get("total_truth_entries") is not None and s.get("n_truth") is not None and int(s["total_truth_entries"]) != int(s["n_truth"]):
                problems.append(f"n_truth {s['n_truth']} != Meta Total_Truth_Entries {s['total_truth_entries']}")
            if pot_check and f.get("pot_used_harvest") and abs(float(s["pot_used"]) - float(f["pot_used_harvest"])) > 1e-6 * abs(float(f["pot_used_harvest"])):
                problems.append(f"pot_used {s['pot_used']} != harvested {f['pot_used_harvest']}")
            f["status"] = "incomplete" if problems else "done"
            f["error"] = "; ".join(problems) if problems else None
            f["pot_used"] = s["pot_used"]
        save_campaign(c)                                   # progress survives an interrupted pass
    save_campaign(c)
    counts: dict[str, int] = {}
    for f in c["files"].values():
        counts[f["status"]] = counts.get(f["status"], 0) + 1
    log(f"campaign {name}: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    return counts


def resubmit(name: str, log=print) -> dict[str, Path]:
    """Worklists of every file not `done`, one per (beam, kind, playlist); returns {worklist: path}."""
    c = load_campaign(name)
    d = campaign_dir(name)
    n = 1 + len([p for p in d.glob("resubmit_*")])
    out = {}
    groups: dict[str, list] = {}
    for tag, f in c["files"].items():
        if f["status"] != "done":
            groups.setdefault(f["worklist"], []).append(f)
    for wl, items in sorted(groups.items()):
        p = d / f"resubmit_{n}_{wl}.txt"
        p.write_text("".join(f["url"] + "\n" for f in sorted(items, key=lambda x: x["line"])))
        for i, f in enumerate(items):
            f["status"] = "planned" if f["status"] in ("failed", "incomplete", "missing") else f["status"]
        out[wl] = p
        log(f"{wl}: {len(items)} files -> {p}")
    save_campaign(c)
    return out


def stage_worklists(name: str, files: list[str] | None = None, log=print) -> dict[str, str]:
    """Upload worklists to `<pnfs_base>/worklists/` so jobs can take them as `-f /pnfs/...` inputs
    (server-side ifdh transfer into $CONDOR_DIR_INPUT). `-f file://` inputs are read by the worker
    itself and fail off-site. Returns {local path: pnfs path} and records it in campaign.json."""
    c = load_campaign(name)
    base = f"{c['pnfs_base']}/worklists"
    todo = files or [w["file"] for w in c["worklists"].values()]
    out = {}
    for f in todo:
        f = str(f)
        dst = f"{base}/{Path(f).name}"
        r = subprocess.run([XRDCP, "-f", "-s", f, f"{DUNE_DOOR}/{xrootd_path(dst)}"], capture_output=True, text=True, timeout=600, env=_token_env())
        if r.returncode != 0:
            raise RuntimeError(f"upload {f} -> {dst}: {r.stderr.strip()[:200]}")
        out[f] = dst
        for w in c["worklists"].values():
            if w["file"] == f:
                w["pnfs_worklist"] = dst
        log(f"staged {Path(f).name} -> {dst}")
    c.setdefault("staged_worklists", {}).update(out)
    save_campaign(c)
    return out


def product_names(tag: str, kind: str, archive: bool = True) -> list[str]:
    """Files a job leaves per AnaTuple: sidecar + reco table (+ MC skims; + the full truth archives when `archive`)."""
    names = [f"manifest_{tag}.json", f"reco_{tag}.npz"]
    if kind == "mc":
        names += [f"truth_{tag}_skim.npz", f"reco_{tag}_truthcols_skim.npz"]
        if archive:
            names += [f"truth_{tag}.npz", f"reco_{tag}_truthcols.npz"]
    return names


def harvest(name: str, products_dir: Path, playlists=None, workers: int = 4, archive: bool = True, log=print) -> dict:
    """Copy every `done` file's products from PNFS into <products_dir>/<beam>/<playlist>/files/<tag>/.

    `archive=False` leaves the full truth tables (`truth_<tag>.npz`, `reco_<tag>_truthcols.npz`, ~260 MB
    per MC file) on PNFS and harvests only what the merge needs (sidecar, reco table, skims, ~65 MB).
    Copies are verified (non-empty; the sidecar must parse and name the tag) and retried; a sidecar that
    stays empty is taken from the campaign's validated `sidecars/` cache.
    """
    from ..io import cheap_fingerprint
    c = load_campaign(name)
    side_dir = campaign_dir(name) / "sidecars"
    todo, n_skip = [], 0
    for tag, f in c["files"].items():
        if f["status"] != "done" or (playlists and f["playlist"] not in playlists):
            continue
        dst = Path(products_dir) / f["beam"] / f["playlist"] / "files" / tag
        if (dst / "harvested.json").exists():
            n_skip += 1; f["harvested"] = str(dst); continue
        todo.append((tag, f, dst))

    def copy_checked(src: str, dst: Path, tag: str):
        for attempt in range(3):
            pnfs_copy(src, dst)
            if dst.stat().st_size > 0:
                if dst.suffix != ".json":
                    return
                try:
                    if json.loads(dst.read_text()).get("tag") == tag:
                        return
                except Exception:  # noqa: BLE001
                    pass
            time.sleep(2 * (attempt + 1))
        if dst.suffix == ".json" and (side_dir / dst.name).exists():
            shutil.copyfile(side_dir / dst.name, dst); return
        raise RuntimeError(f"{src}: copy stayed empty/invalid after 3 attempts")

    def one(item):
        tag, f, dst = item
        names = product_names(tag, f["kind"], archive)
        for n in names:
            copy_checked(f"{f['pnfs_dir']}/{n}", dst / n, tag)
        (dst / "harvested.json").write_text(json.dumps({"harvested": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "from": f["pnfs_dir"],
                                                       "archive": archive, "files": [cheap_fingerprint(dst / n) for n in names]}, indent=1))
        return tag, str(dst)

    n_ok, errors = 0, []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for fu in as_completed([ex.submit(one, it) for it in todo]):
            try:
                tag, d = fu.result()
                c["files"][tag]["harvested"] = d; n_ok += 1
                log(f"harvested {tag} -> {d}")
                if n_ok % 25 == 0:
                    save_campaign(c)
            except Exception as e:  # noqa: BLE001
                errors.append(str(e)[:200]); log(f"harvest failed: {str(e)[:120]}")
    save_campaign(c)
    return {"harvested": n_ok, "already": n_skip, "failed": len(errors), "errors": errors[:10]}
