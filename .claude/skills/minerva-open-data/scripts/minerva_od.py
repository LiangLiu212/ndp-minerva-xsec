#!/usr/bin/env python3
"""minerva_od.py — helpers for the MINERvA Open Data Product (MasterAnaDev AnaTuples on XRootD).

Subcommands
  filelists [--playlist X] [--kind data|mc|<special>] [--beam FHC|RHC] [--out DIR] [--refresh]
      Scrape minerva.fnal.gov/getdata (+ the special-sample pages) for the published per-playlist
      file lists, download them into DIR (default: <skill>/filelists, cached), and print a summary or,
      with filters, the matching lists' xrootd URLs.
  fetch (--list FILE | --playlist X --kind K [--beam B]) --dest DIR [--runs 6038,6040] [--first N] [--dry-run]
      xrdcp every selected file into DIR, skipping files already present with the remote size.
  inspect FILE [FILE ...]
      Trees (with key cycles), entries, branch counts, Meta POT, data/MC classification (needs uproot).
  status [--data-dir DIR]
      AnaTuples and ndp-platform caches on disk (default DIR: ndp.yaml data_dir of the platform).
Kinds: data, mc (=StandardMC), Extended2p2h, CCDiffractivePion, NCDiffractivePion, DSCalorimeters,
ElectronNeutrino, NeutrinoElectronElastic (as named on the site; case-insensitive).
"""
import argparse, os, re, subprocess, sys, time, urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLATFORM = HERE.parents[3]                       # skill/scripts -> skill -> skills -> .claude -> platform
SITE = "https://minerva.fnal.gov"
PAGES = ["getdata", "specialsample2p2h", "specialsamplecohdiffpion", "specialsampledsecal",
         "specialsampleelastic", "specialsamplenue"]
XRD_HOST = "fndcadoor.fnal.gov:1095"                # default door; special-sample lists name fndca1.fnal.gov:1095
LIST_RE = re.compile(r'href="(https://minerva\.fnal\.gov/wp-content/uploads/[^"]+\.txt)"')
NAME_RE = re.compile(r"(MediumEnergy|LowEnergy)_(FHC|RHC)_([A-Za-z0-9]+)_Playlist([A-Za-z0-9]+)(?:-\d+)?\.txt$")
KIND_ALIAS = {"data": "Data", "mc": "StandardMC", "standardmc": "StandardMC"}


def http_get(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "minerva_od.py (curl-like)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def norm_kind(k):
    return KIND_ALIAS.get(k.lower(), k) if k else None


def parse_name(url):
    m = NAME_RE.search(url.rsplit("/", 1)[-1])
    if not m:
        return None
    energy, beam, kind, playlist = m.groups()
    return {"energy": energy, "beam": beam, "kind": kind, "playlist": playlist, "url": url,
            "file": url.rsplit("/", 1)[-1]}


def scrape_lists():
    urls = set()
    for p in PAGES:
        try:
            urls.update(LIST_RE.findall(http_get(f"{SITE}/{p}/")))
        except Exception as e:  # noqa: BLE001
            print(f"[warn] could not read {SITE}/{p}/: {e}", file=sys.stderr)
    entries = [e for e in (parse_name(u) for u in sorted(urls)) if e]
    return entries


def cmd_filelists(a):
    out = Path(a.out) if a.out else HERE / "filelists"
    out.mkdir(parents=True, exist_ok=True)
    index = out / "index.tsv"
    if index.exists() and not a.refresh:
        entries = [dict(zip(("energy", "beam", "kind", "playlist", "file", "url"), l.rstrip("\n").split("\t")))
                   for l in index.read_text().splitlines() if l.strip()]
    else:
        entries = scrape_lists()
        index.write_text("".join(f"{e['energy']}\t{e['beam']}\t{e['kind']}\t{e['playlist']}\t{e['file']}\t{e['url']}\n" for e in entries))
        print(f"scraped {len(entries)} file lists from {len(PAGES)} pages -> {index}")
    kind = norm_kind(a.kind)
    sel = [e for e in entries
           if (not a.playlist or e["playlist"].upper() == a.playlist.upper())
           and (not kind or e["kind"].lower() == kind.lower())
           and (not a.beam or e["beam"].upper() == a.beam.upper())]
    if not (a.playlist or a.kind or a.beam):
        by = {}
        for e in entries:
            by.setdefault((e["energy"], e["beam"], e["kind"]), []).append(e["playlist"])
        for (en, b, k), pls in sorted(by.items()):
            print(f"{en}_{b:3s} {k:24s} {len(pls):2d} playlists: {' '.join(sorted(pls))}")
        print(f"{len(entries)} lists; use --playlist/--kind/--beam to print one list's URLs (downloaded into {out})")
        return 0
    total = 0
    for e in sel:
        local = out / e["file"]
        if not local.exists() or a.refresh:
            local.write_text(http_get(e["url"]))
        lines = [l.strip() for l in local.read_text().splitlines() if l.strip().startswith("root://")]
        total += len(lines)
        print(f"# {e['file']}: {len(lines)} files")
        for l in lines:
            print(l)
    if not sel:
        print("no file list matches the filters", file=sys.stderr); return 1
    print(f"# {total} files in {len(sel)} list(s)", file=sys.stderr)
    return 0


def xrd_size(url):
    m = re.match(r"root://([^/]+)/+(.*)", url)          # lists use fndcadoor.fnal.gov:1095 or fndca1.fnal.gov:1095
    if not m:
        return None
    host, path = m.groups()
    try:
        o = subprocess.run(["xrdfs", host, "stat", "/" + path], capture_output=True, text=True, timeout=60)
        m = re.search(r"Size:\s+(\d+)", o.stdout)
        return int(m.group(1)) if m else None
    except Exception:  # noqa: BLE001
        return None


def cmd_fetch(a):
    if a.list:
        urls = [l.strip() for l in Path(a.list).read_text().splitlines() if l.strip().startswith("root://")]
    else:
        if not (a.playlist and a.kind):
            print("fetch needs --list FILE or --playlist X --kind K", file=sys.stderr); return 2
        ns = argparse.Namespace(playlist=a.playlist, kind=a.kind, beam=a.beam, out=None, refresh=False)
        import io, contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = cmd_filelists(ns)
        if rc:
            return rc
        urls = [l for l in buf.getvalue().splitlines() if l.startswith("root://")]
    if a.runs:
        want = {int(r) for r in a.runs.split(",")}
        urls = [u for u in urls if int(re.search(r"run(\d{8})", u).group(1)) in want]
    if a.first:
        urls = urls[: a.first]
    dest = Path(a.dest); dest.mkdir(parents=True, exist_ok=True)
    total = 0; todo = []
    for u in urls:
        size = xrd_size(u); local = dest / u.rsplit("/", 1)[-1]
        if local.exists() and size is not None and local.stat().st_size == size:
            print(f"skip (present, {size} B): {local.name}"); continue
        todo.append((u, local, size)); total += size or 0
    print(f"{len(todo)} file(s) to fetch, {total/1e9:.2f} GB -> {dest}")
    if a.dry_run:
        for u, local, size in todo:
            print(f"  {size or '?':>12} {u}")
        return 0
    fail = 0
    for u, local, size in todo:
        t0 = time.time()
        r = subprocess.run(["xrdcp", "-f", u, str(local)], capture_output=True, text=True)
        if r.returncode != 0 or (size is not None and local.exists() and local.stat().st_size != size):
            fail += 1; print(f"FAILED {local.name}: {r.stderr.strip()[-200:]}")
        else:
            dt = time.time() - t0
            print(f"ok {local.name}: {local.stat().st_size/1e6:.1f} MB in {dt:.0f}s ({local.stat().st_size/1e6/max(dt,1e-9):.0f} MB/s)")
    return 1 if fail else 0


def cmd_inspect(a):
    try:
        import uproot
    except ImportError:
        print("uproot not importable: run with /opt/pixi/.pixi/envs/default/bin/python3 or inside `pixi run` in ndp-platform", file=sys.stderr)
        return 2
    rc = 0
    for f in a.files:
        try:
            u = uproot.open(f)
        except Exception as e:  # noqa: BLE001
            print(f"== {f}: cannot open ({e})"); rc = 1; continue
        print(f"== {f}  ({os.path.getsize(f)/1e9:.3f} GB)" if os.path.exists(f) else f"== {f}")
        names = {}
        for k in u.keys(recursive=False):
            n, _, cyc = k.partition(";")
            names.setdefault(n, []).append(int(cyc) if cyc else 0)
        for n, cycs in names.items():
            o = u[n]
            if hasattr(o, "num_entries"):
                extra = f"  (key cycles {sorted(cycs)}: highest is read)" if len(cycs) > 1 else ""
                print(f"  {n:14s} {o.num_entries:>10,} entries, {len(o.keys()):>5} branches{extra}")
            else:
                print(f"  {n:14s} {type(o).__name__}")
        kind = "MC" if "Truth" in names else "data"
        if "Meta" in names:
            m = u["Meta"].arrays(library="np")
            pot = {k: float(v.sum()) for k, v in m.items() if k.startswith("POT")}
            print(f"  {kind}: " + ", ".join(f"{k}={v:.6g}" for k, v in pot.items()) +
                  "".join(f", {k}={int(v.sum())}" for k, v in m.items() if k.startswith("Total_")))
        else:
            print(f"  {kind}: no Meta tree")
    return rc


def platform_data_dir():
    y = PLATFORM / "ndp.yaml"
    if y.exists():
        m = re.search(r"^\s*data_dir:\s*(\S+)", y.read_text(), re.M)
        if m:
            return Path(os.path.expandvars(m.group(1).strip("'\"")))
    return None


def cmd_status(a):
    d = Path(a.data_dir) if a.data_dir else (Path(os.environ["NDP_DATA_DIR"]) if os.environ.get("NDP_DATA_DIR") else platform_data_dir())
    if not d or not d.exists():
        print(f"data dir not found: {d}"); return 1
    print(f"== {d}")
    tuples = sorted(d.rglob("MasterAnaDev_*_AnaTuple_run*.root"))
    for t in tuples:
        print(f"  {t.relative_to(d)}  {t.stat().st_size/1e9:8.3f} GB  {time.strftime('%Y-%m-%d', time.localtime(t.stat().st_mtime))}")
    print(f"  {len(tuples)} AnaTuple(s): {sum(1 for t in tuples if '_data_' in t.name)} data, {sum(1 for t in tuples if '_mc_' in t.name)} MC")
    cache = d / "cache"
    if cache.exists():
        for c in sorted(cache.glob("*.npz")):
            print(f"  cache/{c.name}  {c.stat().st_size/1e6:8.1f} MB  {time.strftime('%Y-%m-%d', time.localtime(c.stat().st_mtime))}")
    if (PLATFORM / "ndp").exists():
        print(f"  platform: cd {PLATFORM} && python3 -m ndp data status   (cache freshness per channel)")
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("filelists"); s.add_argument("--playlist"); s.add_argument("--kind"); s.add_argument("--beam")
    s.add_argument("--out"); s.add_argument("--refresh", action="store_true"); s.set_defaults(fn=cmd_filelists)
    s = sub.add_parser("fetch"); s.add_argument("--list"); s.add_argument("--playlist"); s.add_argument("--kind"); s.add_argument("--beam")
    s.add_argument("--dest", required=True); s.add_argument("--runs"); s.add_argument("--first", type=int)
    s.add_argument("--dry-run", action="store_true"); s.set_defaults(fn=cmd_fetch)
    s = sub.add_parser("inspect"); s.add_argument("files", nargs="+"); s.set_defaults(fn=cmd_inspect)
    s = sub.add_parser("status"); s.add_argument("--data-dir"); s.set_defaults(fn=cmd_status)
    a = p.parse_args()
    sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
