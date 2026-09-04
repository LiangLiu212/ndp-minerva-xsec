"""Stage-1 input summary for the CC-inclusive cross-section pipeline.

Streams the AnaTuples declared in a dataset spec (xrootd, no local copies),
and produces a machine-readable summary.json + human-readable summary.md with,
per file: identity, published run period, data-taking timestamps, event counts,
POT, fingerprint checks against the spec, role-specific extras (data/MC), and
a branch-contract check against the analysis branch catalog.

Role-specific "extras" blocks come from the SUMMARIZERS registry — adding more
info for data or MC later means registering a function, not restructuring.
"""
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot
from XRootD import client as xrd_client

import pdg
from runlog_tools import (RunLog, add_label, args_to_inputs, default_outdir,
                          make_parser)

RTOL_FLOAT = 1e-6  # fingerprint tolerance for floats (spec stores 7 sig figs)


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def xrd_stat(url):
    """Remote size + mtime via the xrootd protocol (no data transfer)."""
    m = re.match(r"(root://[^/]+)(/.*)", url)
    if not m:
        return {"error": f"unparseable xrootd url: {url}"}
    fs = xrd_client.FileSystem(m.group(1))
    status, info = fs.stat(m.group(2))
    if not status.ok or info is None:
        return {"error": status.message.strip()}
    return {"size_bytes": info.size, "mtime": info.modtimestr}


def utc_iso(epoch_sec):
    return datetime.fromtimestamp(int(epoch_sec), tz=timezone.utc).isoformat()


def check(name, ok, detail):
    return {"check": name, "pass": bool(ok), "detail": detail}


def compare_expected(expected, observed):
    """Fingerprint comparison: floats with RTOL_FLOAT, ints exact."""
    checks = []
    for key, exp in expected.items():
        obs = observed.get(key)
        if obs is None:
            checks.append(check(key, False, f"expected {exp}, not observed"))
        elif isinstance(exp, float):
            ok = abs(obs - exp) <= RTOL_FLOAT * abs(exp)
            checks.append(check(key, ok, f"expected {exp:.7g}, streamed {obs:.7g}"))
        elif isinstance(exp, list):
            ok = sorted(obs) == sorted(exp)
            checks.append(check(key, ok, f"expected {exp}, streamed {obs}"))
        else:
            checks.append(check(key, obs == exp, f"expected {exp}, streamed {obs}"))
    return checks


# --------------------------------------------------------------------------
# role-specific extras (open-ended registry: add a function, not structure)
# --------------------------------------------------------------------------
def data_extras(f, entry, ctx):
    tree = f["MasterAnaDev"]
    arrs = tree.arrays(["ev_run", "ev_subrun"], library="np")
    runs = np.unique(arrs["ev_run"])
    extras = {
        "ev_run_values": runs.tolist(),
        "run_matches_filename": bool(len(runs) == 1 and int(runs[0]) == entry["run"]),
        "n_subruns": int(len(np.unique(arrs["ev_subrun"]))),
    }
    t = tree.arrays(["ev_gps_time_sec"], library="np")["ev_gps_time_sec"]
    extras["timestamp_span_hours"] = float((t.max() - t.min()) / 3600.0)
    return extras


def mc_extras(f, entry, ctx):
    truth = f["Truth"]
    reco_n = f["MasterAnaDev"].num_entries
    arrs = truth.arrays(["mc_incoming", "mc_current"], library="np")
    n_signal = int(np.count_nonzero((arrs["mc_incoming"] == ctx["nu_mu_pdg"])
                                    & (arrs["mc_current"] == ctx["cc_current"])))
    genie_fams = sorted({b.split("[")[0] for b in truth.keys()
                         if b.startswith("truth_genie_wgt_")})
    cycles = {k: f[k].num_entries for k in f.keys(cycle=True) if k.startswith("Truth;")}
    return {
        "truth_over_reco": float(truth.num_entries / reco_n),
        "truth_signal_rows": n_signal,
        "truth_signal_fraction": float(n_signal / truth.num_entries),
        "n_truth_genie_wgt_families": len(genie_fams),
        "truth_cycles": cycles,
        "live_truth_cycle_used": int(max(int(k.split(";")[1]) for k in cycles)) if cycles else None,
    }


SUMMARIZERS = {"data": data_extras, "mc": mc_extras}


# --------------------------------------------------------------------------
# branch contract
# --------------------------------------------------------------------------
def contract_for_file(f, role, catalog):
    """Check every catalogued branch against this file. Returns per-role results."""
    reco = f["MasterAnaDev"]
    truth = f["Truth"] if (role == "mc" and "Truth" in f) else None
    meta = f["Meta"]
    results = []
    for br in catalog["branches"]:
        name, roles = br["name"], br["roles"]
        if br.get("tree") == "Meta":
            ok = name in meta
            detail = f"Meta/{name} {'present' if ok else 'MISSING'}"
        elif br.get("family"):
            n = len({b.split("[")[0] for b in reco.keys() if b.startswith(name[:-1])})
            if role == "data":
                ok, detail = (n == 0), f"{n} family branches in data reco tree (expect 0)"
            else:
                want = br.get("expected_family_count")
                nt = len({b.split("[")[0] for b in truth.keys() if b.startswith(name[:-1])})
                ok = (n == want and nt == want)
                detail = f"{n} reco / {nt} truth family branches (expect {want})"
        else:
            should_reco = br["in_data"] if role == "data" else br["in_mc_reco"]
            in_reco = name in reco
            ok = (in_reco == should_reco)
            detail = (f"reco tree: {'present' if in_reco else 'absent'}, "
                      f"expected {'present' if should_reco else 'absent'}")
            if in_reco and should_reco:
                tname = reco[name].typename
                if tname != br["dtype"]:
                    ok, detail = False, f"dtype {tname} != catalog {br['dtype']}"
            if role == "mc" and br["in_truth"]:
                in_t = truth is not None and name in truth
                ok = ok and in_t
                detail += f"; truth tree: {'present' if in_t else 'MISSING'}"
        results.append({"branch": name, "roles": roles, "pass": bool(ok), "detail": detail})
    return results


# --------------------------------------------------------------------------
# per-file summary
# --------------------------------------------------------------------------
def summarize_file(entry, playlist, period, catalog, ctx):
    url = entry["url"]
    rec = {"role": entry["role"], "run": entry["run"], "playlist": playlist, "url": url}
    rec["remote"] = xrd_stat(url)

    rec["run_period"] = {
        "window": [period["t_start"], period["t_end"]],
        "run_range": [period["run_first"], period["run_last"]],
        "water_target": period["water_target"], "helium_target": period["helium_target"],
        "beam_config": period["beam_config"],
        # the published run range describes DATA run numbering; MC files use the
        # simulated 11xxxx block, so the check does not apply there
        "run_in_range": (bool(period["run_first"] <= entry["run"] <= period["run_last"])
                         if entry["role"] == "data" else "n/a (MC run numbering)"),
    }

    with uproot.open(url) as f:
        trees = sorted({k.split(";")[0] for k in f.keys(cycle=True)})
        meta = f["Meta"].arrays(["POT_Used", "POT_Total"], library="np")
        observed = {
            "trees": trees,
            "meta_entries": int(len(meta["POT_Used"])),
            "pot_used": float(meta["POT_Used"].sum()),
            "pot_total": float(meta["POT_Total"].sum()),
            "reco_entries": int(f["MasterAnaDev"].num_entries),
        }
        if "Truth" in trees:
            observed["truth_entries"] = int(f["Truth"].num_entries)

        if entry["role"] == "data":
            t = f["MasterAnaDev"].arrays(["ev_gps_time_sec"], library="np")["ev_gps_time_sec"]
            t0, t1 = utc_iso(t.min()), utc_iso(t.max())
            w0 = datetime.fromisoformat(period["t_start"])
            w1 = datetime.fromisoformat(period["t_end"])
            inside = (w0 <= datetime.fromisoformat(t0)) and (datetime.fromisoformat(t1) <= w1)
            rec["timestamps"] = {"first_event_utc": t0, "last_event_utc": t1,
                                 "inside_run_period_window": bool(inside)}
        else:
            rec["timestamps"] = {"note": "simulated (MC timestamps are not physical)"}

        rec["events"] = {k: observed[k] for k in ("reco_entries", "truth_entries") if k in observed}
        rec["pot"] = {k: observed[k] for k in ("pot_used", "pot_total", "meta_entries")}
        rec["fingerprints"] = compare_expected(entry["expected"], observed)
        rec["extras"] = SUMMARIZERS[entry["role"]](f, entry, ctx)
        rec["branch_contract"] = contract_for_file(f, entry["role"], catalog)

    rec["all_fingerprints_pass"] = all(c["pass"] for c in rec["fingerprints"])
    rec["branch_contract_pass"] = all(c["pass"] for c in rec["branch_contract"])
    return rec


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------
def render_md(summary):
    lines = [f"# Input summary — {summary['dataset']}", ""]
    lines.append(f"Generated {summary['generated_utc']} (streaming-only; no local copies).")
    lines.append("")
    lines.append("| role | run | reco events | truth events | POT_Used | first event (UTC) | last event (UTC) | fingerprints | contract |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for f in summary["files"]:
        ts = f["timestamps"]
        lines.append("| {} | {} | {} | {} | {:.6e} | {} | {} | {} | {} |".format(
            f["role"], f["run"], f["events"]["reco_entries"],
            f["events"].get("truth_entries", "—"), f["pot"]["pot_used"],
            ts.get("first_event_utc", "simulated"), ts.get("last_event_utc", "—"),
            "PASS" if f["all_fingerprints_pass"] else "FAIL",
            "PASS" if f["branch_contract_pass"] else "FAIL"))
    lines.append("")
    d = summary["dataset_level"]
    lines.append("## Dataset level")
    lines.append("")
    for k, v in d.items():
        lines.append(f"- **{k}**: {v}")
    lines.append("")
    for f in summary["files"]:
        lines.append(f"## {f['role']} run {f['run']}")
        lines.append("")
        rp = f["run_period"]
        lines.append(f"- run period ({f['playlist']}): runs {rp['run_range'][0]}–{rp['run_range'][1]}, "
                     f"{rp['window'][0]} → {rp['window'][1]}, water {rp['water_target']}, "
                     f"helium {rp['helium_target']}, beam {rp['beam_config']} "
                     f"(run in range: {rp['run_in_range']})")
        if "inside_run_period_window" in f["timestamps"]:
            lines.append(f"- event timestamps inside window: {f['timestamps']['inside_run_period_window']}")
        rm = f["remote"]
        if "size_bytes" in rm:
            lines.append(f"- remote: {rm['size_bytes']:,} bytes, mtime {rm['mtime']}")
        lines.append(f"- extras: {json.dumps(f['extras'], default=str)}")
        failed = [c for c in f["fingerprints"] if not c["pass"]]
        lines.append(f"- fingerprints: {len(f['fingerprints'])-len(failed)}/{len(f['fingerprints'])} pass"
                     + ("" if not failed else f" — FAILED: {failed}"))
        bad = [c for c in f["branch_contract"] if not c["pass"]]
        lines.append(f"- branch contract: {len(f['branch_contract'])-len(bad)}/{len(f['branch_contract'])} pass"
                     + ("" if not bad else f" — FAILED: {bad}"))
        lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------
def main():
    parser = make_parser("Stage-1 input summary: stream AnaTuples, report run period, "
                         "timestamps, events, POT, fingerprints, branch contract.")
    parser.add_argument("--spec", required=True, help="dataset spec JSON (config/datasets/...)")
    parser.add_argument("--catalog", default="config/branches.json",
                        help="analysis branch catalog JSON")
    parser.add_argument("--run-periods", default="config/run_periods.json",
                        help="vendored run-period table JSON")
    parser.add_argument("--outdir", default=None,
                        help="output dir (default: results/<UTC>__summarize_inputs)")
    add_label(parser)
    args = parser.parse_args()

    spec = json.loads(Path(args.spec).read_text())
    catalog = json.loads(Path(args.catalog).read_text())
    periods = json.loads(Path(args.run_periods).read_text())["playlists"]

    outdir = Path(args.outdir) if args.outdir else default_outdir(__file__)
    outdir.mkdir(parents=True, exist_ok=True)

    api = pdg.connect()
    ctx = {
        "nu_mu_pdg": int(api.get_particle_by_name("nu_mu").mcid),
        "cc_current": int(catalog["_meta"]["conventions"]["mc_current"]["CC"]),
    }

    with RunLog(__file__, f"input summary: {spec['_meta']['name']}",
                inputs={**args_to_inputs(args), "dataset": spec["_meta"]["name"],
                        "playlist": spec["playlist"]}) as log:
        playlist = spec["playlist"]
        period = periods[playlist]
        files = [summarize_file(e, playlist, period, catalog, ctx) for e in spec["files"]]

        pot = {f"pot_used_{f['role']}_run{f['run']}": f["pot"]["pot_used"] for f in files}
        pot_data = sum(f["pot"]["pot_used"] for f in files if f["role"] == "data")
        pot_mc = sum(f["pot"]["pot_used"] for f in files if f["role"] == "mc")
        summary = {
            "dataset": spec["_meta"]["name"],
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "playlist": playlist,
            "files": files,
            "dataset_level": {
                "n_files": len(files),
                "runs": sorted(f["run"] for f in files),
                "pot_used_data": pot_data,
                "pot_used_mc": pot_mc,
                "pot_ratio_data_over_mc": (pot_data / pot_mc) if pot_mc else None,
                "reco_events_data": sum(f["events"]["reco_entries"] for f in files if f["role"] == "data"),
                "reco_events_mc": sum(f["events"]["reco_entries"] for f in files if f["role"] == "mc"),
                "truth_events_mc": sum(f["events"].get("truth_entries", 0) for f in files if f["role"] == "mc"),
                "all_fingerprints_pass": all(f["all_fingerprints_pass"] for f in files),
                "branch_contract_pass": all(f["branch_contract_pass"] for f in files),
            },
        }

        json_path = outdir / "summary.json"
        md_path = outdir / "summary.md"
        json_path.write_text(json.dumps(summary, indent=2, default=str))
        md_path.write_text(render_md(summary))

        log.out("outdir", str(outdir))
        log.out("summary_json", str(json_path))
        log.out("summary_md", str(md_path))
        log.out("pot", pot)
        log.out("pot_ratio_data_over_mc", summary["dataset_level"]["pot_ratio_data_over_mc"])
        log.out("all_fingerprints_pass", summary["dataset_level"]["all_fingerprints_pass"])
        log.out("branch_contract_pass", summary["dataset_level"]["branch_contract_pass"])

        print(render_md(summary))
        ok = (summary["dataset_level"]["all_fingerprints_pass"]
              and summary["dataset_level"]["branch_contract_pass"])
        print(f"\n=> overall: {'PASS' if ok else 'FAIL'}  ({json_path})")


if __name__ == "__main__":
    main()
