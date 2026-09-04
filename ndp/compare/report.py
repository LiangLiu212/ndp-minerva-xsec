"""Scorecard JSON + human report for one model run."""
from __future__ import annotations

from pathlib import Path

from ..io import dump_json


def _fmt(x, nd=2):
    return "—" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def write_report(run_dir: Path, ctx: dict) -> Path:
    """ctx keys: model, channel, prediction_summary, unfolded (rows, ranking), folded, figures, warnings."""
    dump_json(ctx, run_dir / "scorecard.json")
    L = []
    m, ch = ctx["model"], ctx["channel"]
    L.append(f"# NDP model test — {m['name']} on {ch['name']}\n")
    L.append(f"*Model kind:* `{m['kind']}`  ·  *Channel:* {ch['description'].strip()}\n")
    if m.get("description"):
        L.append(f"> {m['description'].strip()}\n")
    if ctx.get("prediction_summary"):
        ps = ctx["prediction_summary"]
        L.append("## Model sample\n")
        for k, v in ps.items():
            L.append(f"- **{k}**: {v}")
        L.append("")
    if ctx.get("unfolded"):
        u = ctx["unfolded"]
        L.append("## Unfolded-space comparison (published d²σ + covariance)\n")
        L.append("| prediction | denominator | N cells | χ²_total/ndf | χ²_shape/ndf | α_shape | norm offset % |")
        L.append("|---|---|---|---|---|---|---|")
        for r in u["rows"]:
            L.append(f"| **{r['label']}** | {r['denominator']} | {r['n_populated']} | {_fmt(r['chi2_total_per_ndf'])} | "
                     f"{_fmt(r['chi2_shape_per_ndf'])} | {_fmt(r['alpha_shape'], 3)} | {_fmt(r['norm_offset_pct'], 1)} |")
        for r in u.get("ranking", [])[:20]:
            L.append(f"| {r['label']} | {r['denominator']} | {r['n_populated']} | {_fmt(r['chi2_total_per_ndf'])} | "
                     f"{_fmt(r['chi2_shape_per_ndf'])} | {_fmt(r['alpha_shape'], 3)} | {_fmt(r['norm_offset_pct'], 1)} |")
        L.append("")
        L.append("Rows in bold are this model; the others are the generator curves shipped with the release, scored "
                 "identically for context. Read χ²_shape together with α (a profiled normalisation far from 1 means the "
                 "shape χ² was evaluated at a rescaled model). Total and shape χ² are reported side by side, not combined.\n")
        if u.get("normalisation"):
            L.append(f"Normalisation: {u['normalisation']}\n")
    if ctx.get("folded"):
        f = ctx["folded"]
        g = f["gof"]; t = f["totals"]
        L.append("## Folded-space comparison (model → detector surrogate → reconstructed data)\n")
        L.append(f"- Surrogate: `{f['surrogate']}` ({f.get('surrogate_kind')})")
        L.append(f"- Data: {f['n_data_selected']} selected candidates, {t['data']:.0f} in the grid, POT {f['pot_data']:.4g}")
        L.append(f"- Prediction: {t['pred']:.1f} events ({t['pred_signal']:.1f} signal + {t['bkg']:.1f} background); data/pred = {t['ratio_data_over_pred']:.3f}")
        L.append(f"- Expected true signal: {f['expected']['how']}")
        L.append(f"- Goodness of fit: −2lnL/ndf = {g['minus2lnL']:.1f}/{g['ndf']} ; Pearson χ²/ndf (incl. MC stat) = {g['pearson_chi2']:.1f}/{g['ndf']}")
        if g["data_in_cells_with_zero_prediction"]:
            L.append(f"- {g['data_in_cells_with_zero_prediction']:.0f} data events fall in cells where the prediction is zero (not scored)")
        L.append("")
    if ctx.get("figures"):
        L.append("## Figures\n")
        for fg in ctx["figures"]:
            L.append(f"- `{fg}`")
        L.append("")
    if ctx.get("warnings"):
        L.append("## Caveats\n")
        for w in ctx["warnings"]:
            L.append(f"- {w}")
        L.append("")
    L.append("Provenance: `manifest.json` (inputs, fingerprints, versions, git state) and `scorecard.json` (all numbers).")
    (run_dir / "report.md").write_text("\n".join(L) + "\n")
    return run_dir / "report.md"
