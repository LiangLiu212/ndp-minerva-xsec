"""POT ledger over the ME FHC playlists (streaming Meta-tree reads only).

For every file in every vendored playlist list (config/playlists/), stream
the Meta tree (POT_Used, POT_Total, entry count) and aggregate per playlist
and per role. Compares the per-playlist data sums against the published
getdata-page values (vendored in config/run_periods.json) and the grand
total against the paper's 10.61e20.

Outputs: ledger.json (per-file + per-playlist), ledger.md (tables), RunLog.
"""
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import uproot

from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)

PLAYLISTS = ["1A", "1B", "1C", "1D", "1E", "1F", "1G", "1L", "1M", "1N", "1O", "1P"]


def read_meta(url):
    run = int(re.search(r"run(\d{8})", url).group(1))
    with uproot.open(url) as f:
        meta = f["Meta"].arrays(["POT_Used", "POT_Total"], library="np")
    n = len(meta["POT_Used"])
    assert n >= 1, f"no Meta entries in {url}"
    return {"url": url, "run": run, "meta_entries": n,
            "pot_used": float(meta["POT_Used"].sum()),
            "pot_total": float(meta["POT_Total"].sum())}


def read_with_retry(url, retries=1):
    for attempt in range(retries + 1):
        try:
            return read_meta(url)
        except AssertionError:
            raise
        except Exception as err:
            if attempt == retries:
                return {"url": url, "error": f"{type(err).__name__}: {err}"}
            time.sleep(2.0)


def main():
    parser = make_parser("POT ledger: stream Meta trees of all playlist files, "
                         "aggregate per playlist, compare vs published values.")
    parser.add_argument("--playlists-dir", default="config/playlists")
    parser.add_argument("--run-periods", default="config/run_periods.json")
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--outdir", default=None)
    add_label(parser)
    args = parser.parse_args()

    pdir = Path(args.playlists_dir)
    periods = json.loads(Path(args.run_periods).read_text())["playlists"]
    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    jobs = []  # (playlist, kind, url)
    for pl in PLAYLISTS:
        for kind, tag in (("data", "Data"), ("mc", "StandardMC")):
            lst = pdir / f"MediumEnergy_FHC_{tag}_Playlist{pl}.txt"
            for url in [u.strip() for u in lst.read_text().splitlines() if u.strip()]:
                jobs.append((pl, kind, url))

    with RunLog(__file__, f"POT ledger: {len(jobs)} files, {len(PLAYLISTS)} playlists",
                inputs={**args_to_inputs(args), "n_files": len(jobs)}) as log:
        t0 = time.time()
        results, failures = {}, []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(read_with_retry, url): (pl, kind, url)
                    for pl, kind, url in jobs}
            for i, fut in enumerate(as_completed(futs), 1):
                pl, kind, url = futs[fut]
                r = fut.result()
                if "error" in r:
                    failures.append(r)
                else:
                    results.setdefault((pl, kind), []).append(r)
                if i % 200 == 0 or i == len(jobs):
                    print(f"  [{i}/{len(jobs)}] {len(failures)} failures, "
                          f"{time.time()-t0:.0f}s", flush=True)

        ledger, md = {}, []
        md.append("| playlist | data files | data POT_Used | published (2 s.f.) | ratio | mc files | mc POT_Used | mc/data |")
        md.append("|---|---|---|---|---|---|---|---|")
        tot_d = tot_d_total = tot_m = 0.0
        for pl in PLAYLISTS:
            d = results.get((pl, "data"), [])
            m = results.get((pl, "mc"), [])
            pd = sum(r["pot_used"] for r in d)
            pm = sum(r["pot_used"] for r in m)
            pub = periods[f"minervame{pl}"]["data_pot_e20"]
            tot_d += pd
            tot_d_total += sum(r["pot_total"] for r in d)
            tot_m += pm
            ledger[f"minervame{pl}"] = {
                "data": {"n_files": len(d), "pot_used": pd,
                         "pot_total": sum(r["pot_total"] for r in d),
                         "published_e20": pub},
                "mc": {"n_files": len(m), "pot_used": pm},
                "files": {"data": d, "mc": m},
            }
            md.append(f"| {pl} | {len(d)} | {pd:.4e} | {pub}e20 | "
                      f"{pd/(pub*1e20):.3f} | {len(m)} | {pm:.4e} | "
                      f"{pm/pd:.2f} |")
        md.append(f"| **total** | {sum(len(results.get((p,'data'),[])) for p in PLAYLISTS)} "
                  f"| **{tot_d:.4e}** |  |  | "
                  f"{sum(len(results.get((p,'mc'),[])) for p in PLAYLISTS)} "
                  f"| **{tot_m:.4e}** | {tot_m/tot_d:.2f} |")

        comparison = {
            "total_data_pot_used": tot_d,
            "total_data_pot_total": tot_d_total,
            "total_mc_pot_used": tot_m,
            "paper_pot_e20": 10.61,
            "ratio_to_paper": tot_d / 10.61e20,
            "getdata_page_sum_e20": 11.12,
            "ratio_to_page_sum": tot_d / 11.12e20,
            "n_failures": len(failures),
            "failed_urls": [f["url"] for f in failures],
            "wall_s": round(time.time() - t0, 1),
        }

        (outdir / "ledger.json").write_text(json.dumps(
            {"per_playlist": ledger, "comparison": comparison}, indent=2))
        md_text = "\n".join(md) + "\n\n" + json.dumps(comparison, indent=2)
        (outdir / "ledger.md").write_text(md_text)

        log.out("outdir", str(outdir))
        log.out("ledger_json", str(outdir / "ledger.json"))
        log.out("comparison", comparison)
        print("\n".join(md))
        print(json.dumps(comparison, indent=2))


if __name__ == "__main__":
    main()
