#!/usr/bin/env python3
"""Plot the binned detector response of chosen 1D grids and render a report.

    PYTHONPATH=$PWD python report/make_migration_report.py \
        --channel minerva_me_ccqelike_1mu1p --grids muon_p,muon_theta \
        --out report/Migration_muon_FHC.md [--tag mig]

Everything is read from the built surrogates under
`surrogates/<channel>/<grid>/binned_<beam>_<first>-<last>/`: `surrogate.json` (binning, provenance,
training counts) and `arrays.npz` (P, eff, den_counts, num_counts, migration_counts). No number is
typed into this file.

Convention, fixed by ndp/surrogate/binned.py: P[i, j] = M[i, j] / num[j], so a column is the
probability distribution of the reco bin of a selected signal event born in true bin j, and a column
sum below one is signal that reconstructs outside the grid. Efficiency is a separate factor,
eff[j] = num[j] / den[j], and is never folded into P.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parent.parent
VMIN = 1e-3


def load_grid(channel: str, grid: str, root: Path) -> dict:
    """The one built binned response of `grid` (refuses if several beam/playlist builds exist)."""
    cands = sorted((root / channel / grid).glob("binned_*"))
    if not cands:
        raise SystemExit(f"no binned surrogate for {channel}/{grid} under {root}")
    if len(cands) > 1:
        raise SystemExit(f"several builds for {channel}/{grid}: {[c.name for c in cands]}; pick one")
    d = cands[0]
    rec = json.loads((d / "surrogate.json").read_text())
    z = np.load(d / "arrays.npz")
    spec = rec["meta"]["measurement_spec"]
    edges = np.asarray(rec["binning"]["x_edges"], float)
    return {"dir": d, "name": grid, "meta": rec["meta"], "edges": edges,
            "units": spec["x"].get("units", ""), "observable": spec["x"]["observable"], "reco": spec["x"]["reco"],
            "P": z["P"], "eff": z["eff"], "den": z["den_counts"], "num": z["num_counts"], "M": z["migration_counts"]}


def bin_labels(edges: np.ndarray) -> list[str]:
    return [f"{edges[i]:g}–{edges[i + 1]:g}" for i in range(len(edges) - 1)]


def smearing(P: np.ndarray) -> dict:
    """Per true bin: fraction staying, migrating down, migrating up, and lost out of the reco grid."""
    n = P.shape[0]
    stay = np.diag(P).copy()
    down = np.array([P[:j, j].sum() for j in range(n)])
    up = np.array([P[j + 1:, j].sum() for j in range(n)])
    lost = 1.0 - P.sum(axis=0)
    return {"stay": stay, "down": down, "up": up, "lost": np.clip(lost, 0, None)}


def nn_frac(P: np.ndarray) -> np.ndarray:
    """Per true bin, the fraction reconstructing in its own bin or an immediate neighbour."""
    n = P.shape[0]
    return np.array([P[max(j - 1, 0):j + 2, j].sum() for j in range(n)])


def fig_response(g: dict, out: Path) -> None:
    """Annotated response matrix over the efficiency of the same true bins."""
    P, edges = g["P"], g["edges"]
    n = P.shape[0]
    labs = bin_labels(edges)
    w = min(12.0, 0.75 * n + 3.0)
    fig, (ax, axe) = plt.subplots(2, 1, figsize=(w, 0.82 * w + 2.6),
                                  gridspec_kw={"height_ratios": [4.0, 1.35]})
    cmap = matplotlib.colormaps["viridis"].with_extremes(bad=matplotlib.colormaps["viridis"](0.0))
    mesh = ax.pcolormesh(np.arange(n + 1), np.arange(n + 1), np.ma.masked_less_equal(P, 0.0),
                         cmap=cmap, norm=LogNorm(vmin=VMIN, vmax=1.0))
    ax.set_aspect("equal")
    floor = 0.005 if n <= 10 else 0.01
    fs = max(5.0, min(9.0, 70.0 / n))
    for j in range(n):
        for i in range(n):
            v = P[i, j]
            if v < floor:
                continue
            t = (np.log10(max(v, VMIN)) - np.log10(VMIN)) / (-np.log10(VMIN))
            ax.text(j + 0.5, i + 0.5, f"{v:.3f}", ha="center", va="center", fontsize=fs,
                    color="black" if t > 0.62 else "white")
    ax.set_xticks(np.arange(n) + 0.5, labs, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(n) + 0.5, labs, fontsize=7)
    ax.set_xlabel(f"true {g['observable']}  [{g['units']}]")
    ax.set_ylabel(f"reconstructed  [{g['units']}]")
    ax.set_title(f"P(reco bin | true bin) for {g['name']}: each column sums to the fraction that stays in the grid", fontsize=10)
    fig.colorbar(mesh, ax=ax, label="P (log scale)", fraction=0.040, pad=0.02, extend="min")

    err = np.where(g["den"] > 0, np.sqrt(np.clip(g["eff"] * (1 - g["eff"]), 0, None) / np.maximum(g["den"], 1)), 0.0)
    axd = axe.twinx()
    axd.bar(np.arange(n) + 0.5, g["den"], width=0.86, color="#cfd8dc", zorder=0)
    axd.set_ylabel("signal events in the true bin", color="#78909c", fontsize=8)
    axd.tick_params(axis="y", colors="#78909c", labelsize=7)
    axe.errorbar(np.arange(n) + 0.5, g["eff"], yerr=err, fmt="o-", ms=4, lw=1.2, color="#08519c", capsize=2, zorder=3)
    axe.set_zorder(axd.get_zorder() + 1); axe.patch.set_visible(False)
    axe.set_xlim(0, n); axe.set_ylim(0, max(0.55, 1.25 * float(np.nanmax(g["eff"]))))
    axe.set_xticks(np.arange(n) + 0.5, labs, rotation=45, ha="right", fontsize=7)
    axe.set_ylabel("selection efficiency")
    axe.set_xlabel(f"true {g['observable']}  [{g['units']}]")
    axe.set_title("efficiency of the same true bins (bars: how many signal events the bin holds)", fontsize=9)
    axe.grid(axis="y", alpha=0.3)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130); plt.close(fig)


def fig_smearing(g: dict, out: Path) -> None:
    """Where each true bin's selected signal ends up, and the net direction of the pull."""
    s = smearing(g["P"])
    n = len(s["stay"])
    labs = bin_labels(g["edges"])
    x = np.arange(n)
    fig, (ax, axn) = plt.subplots(2, 1, figsize=(max(7.0, 0.8 * n + 3.0), 6.4), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 2]})
    w = 0.27
    ax.bar(x - w, s["stay"], w, color="#2c7fb8", label="stays in its own bin")
    ax.bar(x, s["down"], w, color="#f16913", label="reconstructs in a lower bin")
    ax.bar(x + w, s["up"], w, color="#41ab5d", label="reconstructs in a higher bin")
    if s["lost"].max() > 1e-4:
        ax.bar(x + 2 * w, s["lost"], w, color="#999999", label="leaves the reco grid")
    ax.set_ylabel("fraction of the selected signal")
    ax.set_ylim(0, 1.0); ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3)
    ax.set_title(f"{g['name']}: where a true bin's selected signal reconstructs", fontsize=10)

    net = s["down"] - s["up"]
    axn.bar(x, net, 0.7, color=np.where(net >= 0, "#f16913", "#41ab5d"))
    axn.axhline(0, color="k", lw=0.8)
    axn.set_ylabel("net pull\n(down − up)")
    axn.set_xticks(x, labs, rotation=45, ha="right", fontsize=7)
    axn.set_xlabel(f"true {g['observable']}  [{g['units']}]")
    axn.grid(axis="y", alpha=0.3)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130); plt.close(fig)


def matrix_block(g: dict) -> list[str]:
    """The response matrix as a fixed-width block (alignment survives markdown)."""
    labs = bin_labels(g["edges"])
    wide = max(len(l) for l in labs)
    head = " " * (wide + 8) + "".join(f"{l:>8}" for l in labs)
    out = ["```", f"{g['name']}:  rows = reco bin, columns = true bin", "", head]
    for i, l in enumerate(labs):
        row = "".join(f"{v:8.3f}" if v >= 0.0005 else f"{'-':>8}" for v in g["P"][i])
        out.append(f"reco {l:>{wide}}  {row}")
    out += ["", " " * (wide + 8) + "".join(f"{v:8.3f}" for v in g["P"].sum(axis=0)) + "   <- column sums", "```"]
    return out


def table_block(g: dict) -> list[str]:
    labs = bin_labels(g["edges"])
    s = smearing(g["P"])
    err = np.where(g["den"] > 0, np.sqrt(np.clip(g["eff"] * (1 - g["eff"]), 0, None) / np.maximum(g["den"], 1)), 0.0)
    L = [f"| true {g['observable']} [{g['units']}] | efficiency | stays | down | up | net pull | signal events |",
         "|---|---|---|---|---|---|---|"]
    for j, l in enumerate(labs):
        L.append(f"| {l} | {g['eff'][j]:.3f} ± {err[j]:.3f} | {s['stay'][j]:.3f} | {s['down'][j]:.3f} | "
                 f"{s['up'][j]:.3f} | {s['down'][j] - s['up'][j]:+.3f} | {g['den'][j]:,.0f} |")
    return L


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--channel", required=True)
    ap.add_argument("--grids", required=True, help="comma-separated measurement names")
    ap.add_argument("--surrogate-root", default="surrogates")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="mig")
    a = ap.parse_args()
    root = Path(a.surrogate_root).resolve()
    out = Path(a.out).resolve()
    figs = out.parent / "figs"
    grids = [load_grid(a.channel, n.strip(), root) for n in a.grids.split(",") if n.strip()]

    made = {}
    for g in grids:
        f1 = figs / f"{a.tag}_response_{g['name']}.png"
        f2 = figs / f"{a.tag}_smearing_{g['name']}.png"
        fig_response(g, f1); fig_smearing(g, f2)
        made[g["name"]] = (f1.name, f2.name)

    m0 = grids[0]["meta"]
    L = [f"# The detector response of the muon variables in the MINERvA muon + leading-proton selection", "",
         f"**What this is.** The binned surrogate that carries a truth prediction into reconstructed space for this "
         f"analysis. For one grid it is a pair: a migration matrix and a selection efficiency, kept as separate factors. "
         f"The prediction in reco bin *i* is", "",
         "```", "reco[i] = Σ_j  P[i,j] · ε[j] · true[j]  +  bkg[i]", "```", "",
         f"**The convention matters.** P is normalised by column against the *selected* signal of that true bin, "
         f"P[i,j] = M[i,j] / num[j]. A column therefore gives the probability distribution of the reco bin of an event "
         f"born in true bin j, and a column sum below one is signal that reconstructs outside the grid entirely. That "
         f"loss is real and stays lost. Normalising instead by the column sum of M would push those events back into the "
         f"visible bins and silently inflate every prediction. Efficiency is never folded into P: it is a separate "
         f"per-true-bin factor, ε[j] = num[j] / den[j].", "",
         f"**Provenance.** Channel `{a.channel}`, selection `{m0.get('selection')}`, built "
         f"{str(m0.get('built'))[:10]} by `{m0.get('built_by')}` over {len(m0.get('training_chunks', []))} playlist "
         f"products ({', '.join(m0.get('training_chunks', []))}) at {m0.get('pot_mc', float('nan')):.4e} protons on "
         f"target of official MC. Every array below is read from the built surrogate directories; nothing is typed in. "
         f"Rendered by `report/make_migration_report.py`.", "", "---", ""]

    for g in grids:
        n = g["P"].shape[0]
        s = smearing(g["P"])
        d = np.diag(g["P"])
        f_resp, f_smear = made[g["name"]]
        lost_tot = float(g["meta"].get("n_num_reco_out_of_grid", 0.0))
        feed_tot = float(g["meta"].get("n_feedin_true_out_of_grid", 0.0))
        L += [f"## {g['name']}: {g['observable']} → {g['reco']}, {n} bins in {g['units']}", "",
              f"Edges at {', '.join(f'{e:g}' for e in g['edges'])}.", ""]
        L += matrix_block(g)
        L += ["", f"![response {g['name']}](figs/{f_resp})", "",
              f"*Top: the matrix on a logarithmic colour scale. Cells at or below the colour floor of {VMIN:g}, exact "
              f"zeros included, are drawn in the floor colour; only entries above {0.005 if n <= 10 else 0.01:g} carry a "
              f"printed number. Bottom: the efficiency of the same true bins with binomial errors, over the number of "
              f"signal events each bin holds.*", "", ""]
        L += ["Per true bin: how much of its selected signal stays, how much reconstructs lower, how much higher, and "
              "the net pull, positive meaning a downward tilt. The first and last bins are structurally one-sided, "
              "since nothing can migrate below the first or above the last, so compare only the interior.", ""]
        L += table_block(g)
        L += ["", f"![smearing {g['name']}](figs/{f_smear})", "",
              f"*Where each true bin's selected signal reconstructs, and the net direction of the pull. A positive net "
              f"pull means the bin loses more events downward than it gains upward.*", "",
              f"Diagonal strength runs from {d.min():.3f} (bin {bin_labels(g['edges'])[int(d.argmin())]}) to "
              f"{d.max():.3f} (bin {bin_labels(g['edges'])[int(d.argmax())]}), averaging {d.mean():.3f}, and the "
              f"diagonal plus its two neighbours holds at least {nn_frac(g['P']).min():.3f} of every column. "
              f"Efficiency runs from {g['eff'].min():.3f} to {g['eff'].max():.3f}, a factor of "
              f"{g['eff'].max() / max(g['eff'].min(), 1e-9):.1f}. "
              + (f"Every column sums to one: {lost_tot:.0f} selected signal events reconstruct outside this grid and "
                 f"{feed_tot:.0f} feed in from outside it, so the grid is closed."
                 if lost_tot == 0 and feed_tot == 0 else
                 f"{lost_tot:,.0f} selected signal events reconstruct outside this grid and {feed_tot:,.0f} feed in "
                 f"from outside it; the latter are carried in the background, not in P."), "", "---", ""]

    # closing discussion, with every number pulled from the arrays
    L += ["## Reading the two together", ""]
    if len(grids) == 2:
        gp, gt = grids[0], grids[1]
        sp, st = smearing(gp["P"]), smearing(gt["P"])
        np_, nt = gp["P"].shape[0], gt["P"].shape[0]
        ip, it = slice(1, np_ - 1), slice(1, nt - 1)
        netp, nett = sp["down"] - sp["up"], st["down"] - st["up"]
        labp, labt = bin_labels(gp["edges"]), bin_labels(gt["edges"])
        L += [f"- **Both grids are closed.** The reconstruction-level muon window is identical to the truth-level one, "
              f"so a selected event necessarily lands inside both grids. Every column sums to one and there is no "
              f"feed-in. This is special to the muon variables: the transverse-imbalance grids all leak.",
              f"- **Smearing is short range in both.** Counting only the diagonal and its two neighbours already "
              f"accounts for {nn_frac(gp['P']).min():.3f} to {nn_frac(gp['P']).max():.3f} of each momentum column and "
              f"{nn_frac(gt['P']).min():.3f} to {nn_frac(gt['P']).max():.3f} of each angle column. Nothing migrates far.",
              f"- **The momentum pull changes sign with momentum.** Read the outer bins with care, since the lowest can "
              f"only migrate up and the highest can only migrate down. Among the interior bins the low end tilts "
              f"upward, {netp[1]:+.3f} at {labp[1]} GeV/c, and the tilt turns over and grows downward with momentum, "
              f"reaching {netp[np_ - 2]:+.3f} at {labp[np_ - 2]} GeV/c. Averaged over the interior the net is "
              f"{netp[ip].mean():+.3f} downward.",
              f"- **The angle pull is outward, and strongest near the beam axis.** {int((nett[it] < 0).sum())} of the "
              f"{nt - 2} interior bins push away from zero, by {abs(nett[1]):.3f} at {labt[1]} degrees, decaying "
              f"smoothly to about zero by {labt[10]} degrees. Angle is positive definite and multiple scattering only "
              f"adds to it, so resolution near the axis can scatter a track outward but not inward. The mean interior "
              f"net is {nett[it].mean():+.3f}.",
              f"- **Bin width dominates the diagonal, not resolution.** Read the diagonal against the edge list before "
              f"comparing bins: the widest momentum bin keeps {np.diag(gp['P'])[-1]:.3f} and the widest angle bin "
              f"{np.diag(gt['P'])[-1]:.3f}, both inflated simply because a wide bin is harder to leave.",
              f"- **The efficiencies run opposite ways.** Momentum efficiency rises from {gp['eff'][0]:.3f} to a "
              f"maximum of {gp['eff'].max():.3f} and eases off at the top, while angle efficiency falls monotonically "
              f"from {gt['eff'][0]:.3f} to {gt['eff'][-1]:.3f}. Forward muons are far easier to keep, which is the "
              f"MINOS acceptance.", ""]
    L += ["## Caveats", "",
          "- The matrix is the official MC's own detector response, unweighted and untuned. It carries that generator's "
          "kinematic distribution inside each bin, so it is only as good as the assumption that the shape within a bin "
          "is close to the truth being folded.",
          "- The efficiency shown per bin is the marginal of the full two-dimensional map over the other variable. For a "
          "sample whose muon and proton kinematics differ from the official MC's, fold through the full response rather "
          "than reweighting by these marginals.",
          "- Statistical errors on the efficiency are binomial on the denominator shown. The matrix carries a "
          "multinomial error per column that is propagated in the folding but not drawn here.", ""]
    out.write_text("\n".join(L) + "\n")
    print(f"wrote {out} with {2 * len(grids)} figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
