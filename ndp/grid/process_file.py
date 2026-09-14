"""Process one AnaTuple (local path or root:// URL) exactly as `ndp data cache` does, plus the
channel's truth-level signal cutflow, reco selection cutflow, derived-column skims and a
sidecar manifest — the unit of work of a grid job, and the same entry point for local parity runs.

    python -m ndp.grid.process_file --url root://... --channel minerva_me_ccqelike_1mu1p --kind mc --out DIR

Exit codes: 0 ok, 2 bad arguments, 5 processing failed after retries (a `manifest_<tag>.json`
with status "failed" is still written so the campaign can see why).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

DOORS = ("fndcadoor.fnal.gov:1095", "fndca1.fnal.gov:1095")
BACKOFF_S = (30, 120, 300)


def swap_door(url: str, door: str) -> str:
    """Rewrite the host:port of a root:// URL."""
    if not url.startswith("root://"):
        return url
    rest = url[len("root://"):]
    _, _, path = rest.partition("/")
    return f"root://{door}/{path}"


def process(url: str, channel_name: str, kind: str, out_dir: str | Path, *, attempts: int = 3, timeout: float = 300.0,
            skim: bool = True, entry_stop: int | None = None, extra_meta: dict | None = None, log=print) -> dict:
    from ..adapters.minerva_anatuple import build_cache, cache_tag
    from ..channels import load_channel
    ch = load_channel(channel_name)
    is_mc = kind == "mc"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = cache_tag(url)
    errors = []
    for k in range(attempts):
        door = DOORS[k % len(DOORS)]
        u = swap_door(url, door) if url.startswith("root://") else url
        t0 = time.time()
        try:
            res = build_cache(u, out_dir, is_mc=is_mc, truth=is_mc, reco=True, entry_stop=entry_stop, log=log,
                              channel=ch, skim=(skim and is_mc), sidecar=True, timeout=timeout,
                              extra_meta={"attempt": k + 1, "door": door, "requested_url": url, **(extra_meta or {})})
            res["attempts"] = k + 1
            return res
        except Exception as e:  # noqa: BLE001 — any I/O or decode failure is retried
            err = f"attempt {k + 1} via {door}: {type(e).__name__}: {str(e)[:300]}"
            errors.append(err)
            log(err)
            log(traceback.format_exc(limit=3))
            if k + 1 < attempts:
                time.sleep(BACKOFF_S[min(k, len(BACKOFF_S) - 1)])
    side = {"tag": tag, "source": url, "kind": kind, "status": "failed", "errors": errors, "attempts": attempts,
            "channel": channel_name, "built": time.strftime("%Y-%m-%dT%H:%M:%S%z"), **(extra_meta or {})}
    (out_dir / f"manifest_{tag}.json").write_text(json.dumps(side, indent=2))
    return {"tag": tag, "path": url, "status": "failed", "errors": errors}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--url", required=True, help="root:// URL or local path of the AnaTuple")
    ap.add_argument("--channel", required=True)
    ap.add_argument("--kind", choices=("data", "mc"), required=True)
    ap.add_argument("--out", required=True, help="output directory (the job's scratch)")
    ap.add_argument("--attempts", type=int, default=3)
    ap.add_argument("--timeout", type=float, default=300.0, help="XRootD request timeout [s]")
    ap.add_argument("--no-skim", action="store_true")
    ap.add_argument("--entry-stop", type=int, default=None, help="debug: read only the first N entries")
    a = ap.parse_args(argv)
    extra = {k: os.environ.get(v) for k, v in (("cluster", "CLUSTER"), ("process", "PROCESS"), ("payload", "NDP_PAYLOAD"),
                                                 ("host", "HOSTNAME"), ("ndp_git_sha", "NDP_GIT_SHA"))}
    res = process(a.url, a.channel, a.kind, a.out, attempts=a.attempts, timeout=a.timeout, skim=not a.no_skim,
                  entry_stop=a.entry_stop, extra_meta=extra)
    print(json.dumps({k: v for k, v in res.items() if k not in ("fingerprint",)}, indent=1, default=str))
    return 5 if res.get("status") == "failed" else 0


if __name__ == "__main__":
    sys.exit(main())
