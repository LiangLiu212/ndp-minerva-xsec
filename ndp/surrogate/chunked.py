"""Binned surrogates for several measurements in one pass over the MC, one playlist at a time.

`build.build_surrogates` loads the whole MC (fine for one file); the playlist products of a campaign are
~90 GB of arrays, so this module walks `products.iter_mc_chunks` once, evaluates every requested
measurement on each chunk and accumulates the additive counts of `binned.count_pairs` (den, num,
migration, feed-in, background per reco cell and per background category). The result is bit-for-bit
what `BinnedResponse.fit` returns on the concatenated tables, and the closure is checked from the
accumulated counts: fold(den) + background == the selected reco histogram of the same MC.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..io import timestamp
from .binned import BinnedResponse, add_counts, count_pairs

BKG_CATEGORIES = ("bkg 1 pi+-", "bkg 1 pi0", "bkg multi-pi", "bkg other (no pion)")   # diagnostics._CATEGORIES background part


def _tag(labels: list[str]) -> str:
    """`FHC_1A-1P` for playlist products, the comma-joined tags for legacy caches."""
    if labels and all("/" in lab for lab in labels):
        beams = sorted({lab.split("/", 1)[0] for lab in labels})
        pls = [lab.split("/", 1)[1] for lab in labels]
        return f"{'+'.join(beams)}_{pls[0]}-{pls[-1]}" if len(pls) > 1 else f"{beams[0]}_{pls[0]}"
    return ",".join(labels)


def accumulate_counts(channel, measurements: list, cfg, log=print) -> dict:
    """One pass over the MC chunks -> per measurement the summed counts, the selected reco histogram and the
    truth cells; plus POT, labels, sources and the global training counts."""
    from ..channels.selections import select
    from ..diagnostics import mc_categories
    from ..products import iter_mc_chunks
    counts = {m.name: None for m in measurements}
    reco_sel = {m.name: np.zeros(m.binning.n_cells) for m in measurements}
    truth_cells = {m.name: np.zeros(m.binning.n_cells) for m in measurements}
    tot = {"n_den": 0, "n_num": 0, "n_bkg": 0, "n_passed": 0}
    pot_mc = 0.0; labels, sources = [], []
    for label, rm, rt, truth, pot, srcs in iter_mc_chunks(cfg, channel):
        sig_den = channel.is_signal(truth) & channel.in_phase_space(truth)
        passed = select(channel, rm)
        sig_rt = channel.is_signal(rt) & channel.in_phase_space(rt)
        sig_num = passed & sig_rt
        bkg = passed & ~sig_rt
        cats = np.asarray(mc_categories(channel, rt))[bkg]
        tot["n_den"] += int(sig_den.sum()); tot["n_num"] += int(sig_num.sum()); tot["n_bkg"] += int(bkg.sum()); tot["n_passed"] += int(passed.sum())
        # truth observables per measurement on the Truth skim, then on the reco-side truth, then reco (keeps the TKI caches warm)
        td = {m.name: tuple(a[sig_den] for a in m.truth_observables(channel, truth)) for m in measurements}
        tn = {m.name: tuple(a[sig_num] for a in m.truth_observables(channel, rt)) for m in measurements}
        for m in measurements:
            xr, yr = m.reco_observables(rm, params=channel.observable_params)
            c = count_pairs(m.binning, x_true_den=td[m.name][0], y_true_den=td[m.name][1],
                            x_true_num=tn[m.name][0], y_true_num=tn[m.name][1], x_reco_num=xr[sig_num], y_reco_num=yr[sig_num],
                            x_reco_bkg=xr[bkg], y_reco_bkg=yr[bkg], bkg_category=cats, category_names=BKG_CATEGORIES)
            counts[m.name] = add_counts(counts[m.name], c)
            reco_sel[m.name] += m.binning.histogram(xr[passed], yr[passed])[0]
            truth_cells[m.name] += m.binning.histogram(td[m.name][0], td[m.name][1], truth["weight"][sig_den])[0]
        pot_mc += float(pot); labels.append(label); sources += list(srcs)
        log(f"{label}: den {int(sig_den.sum())}, selected {int(passed.sum())} (signal {int(sig_num.sum())}, bkg {int(bkg.sum())}), POT {pot:.4g}")
        del rm, rt, truth, td, tn
    return {"counts": counts, "reco_selected": reco_sel, "truth_cells": truth_cells, "pot_mc": pot_mc, "labels": labels,
            "sources": sources, "training_counts": tot}


def build_surrogates_chunked(channel, measurements: list, cfg, out_root: Path | None = None, log=print) -> dict:
    """Build and save one binned response per measurement from a single pass over the MC.

    Returns {measurement name: {"path", "diagnostics", "closure"}}. Closure is exact by construction
    (fold(den) + bkg_per_pot x pot_mc == selected reco histogram); a non-zero deviation is a bug.
    """
    from ..compare import plots
    from ..products import sources_fingerprints
    acc = accumulate_counts(channel, measurements, cfg, log=log)
    tag = _tag(acc["labels"])
    fps = sources_fingerprints(cfg, channel)
    out = {}
    for m in measurements:
        c = acc["counts"][m.name]
        meta = {"channel": channel.name, "measurement": m.name, "measurement_spec": m.to_dict(), "built": timestamp(),
                "training_mc": fps, "training_sources": acc["sources"], "training_chunks": acc["labels"], "pot_mc": acc["pot_mc"],
                "generator": "MINERvA official MC (playlist products)" if "/" in acc["labels"][0] else "MINERvA official MC",
                "selection": channel.selection.get("name"), "phase_space": channel.phase_space, "signal": channel.signal,
                "training_counts": acc["training_counts"], "built_by": "ndp.surrogate.chunked.build_surrogates_chunked"}
        sur = BinnedResponse.from_counts(m.binning, c, pot_mc=acc["pot_mc"], meta=meta)
        root = Path(out_root) if out_root else m.surrogate_root(cfg)
        d = root / f"binned_{tag}"
        sur.save(d)
        plots.response_figure(sur, d / "response.png", f"binned response, MC {tag}, {m.name}")
        pred = sur.fold(acc["truth_cells"][m.name]) + sur.background(acc["pot_mc"])
        h = acc["reco_selected"][m.name]
        dev = pred - h
        clo = {"n_pred": float(pred.sum()), "n_reco": float(h.sum()), "max_abs_dev": float(np.abs(dev).max()),
               "exact": bool(np.allclose(pred, h, atol=1e-6)), "compared_with": "all selected candidates",
               "ratio_pred_over_reco": float(pred.sum() / h.sum()) if h.sum() else None}
        out[m.name] = {"kind": "binned", "path": str(d), "diagnostics": sur.diagnostics(), "closure": clo}
        log(f"{m.name}: saved {d}; closure max|pred-reco| = {clo['max_abs_dev']:.3g} on {clo['n_reco']:.0f} selected in-grid events"
            + ("" if clo["exact"] else "  ** NOT EXACT **"))
    return out
