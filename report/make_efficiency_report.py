#!/usr/bin/env python3
"""Publish an `ndp efficiency run` directory as a report: its report.md with the figures copied next to it.

    python report/make_efficiency_report.py --run runs/<date>_efficiency_<channel> --out report/Efficiency_1mu1p_FHC.md [--tag eff]

Every number comes from the run directory (`eff_<grid>.json`, `background_<grid>.json`, `ansatz_closure.json`,
`summary.json`); re-render the run's own report first with `python -m ndp efficiency report --run <dir>` if its
tables changed.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

HEAD = """# Selection efficiency maps and background of the MINERvA 1μ1p selection (full ME FHC official MC)

**Status:** rendered {today} from `{run}` (`ndp efficiency run`), which reads the binned responses built by
`ndp surrogate build --channel {channel} --measurement all` over the 12 FHC playlists
(`surrogates/{channel}/<grid>/binned_FHC_1A-1P/`, every closure exact). Per-grid tables live in that run directory:
`eff_<grid>.json` / `.npz` (edges, efficiency, binomial error, denominator, numerator), `background_<grid>.json`
(total, by category, feed-in, at the data POT) and `ansatz_closure.json` (the per-bin closure of the factorised maps).

**How to use this on a truth sample.** For a weight in the bins of one released variable, take that grid's per-bin ε
below (exact for its own binning). For a per-event weight from the two maps,

```
python -m ndp efficiency apply --channel {channel} --run {run} --sample <truth.npz> --out weights.npz
```

writes w = ε_μ(p_μ, cos θ_μ) · ε_p(p_p, cos θ_p) / ⟨ε⟩ per event (0 outside the maps or without a leading proton in
the window). The closure table says how well that factorised weight reproduces the actually selected signal on each
grid: within a few per cent on the muon and proton kinematics, but far off on the transverse-imbalance variables,
where the selection efficiency depends on the muon–proton correlation. For those, fold through the full response
(`python -m ndp run <model> --channel {channel} --measurement all --modes folded`).

---

"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True); ap.add_argument("--out", required=True); ap.add_argument("--tag", default="eff")
    a = ap.parse_args()
    run = Path(a.run).resolve(); out = Path(a.out).resolve(); figs = out.parent / "figs"
    body = (run / "report.md").read_text()
    summary = json.loads((run / "summary.json").read_text())
    figs.mkdir(parents=True, exist_ok=True)
    kept = []
    for src in sorted((run / "figs").glob("*.png")):
        dst = figs / f"{a.tag}_{src.name}"
        shutil.copyfile(src, dst); kept.append(dst.name)
    body = re.sub(r"^# .*\n", "", body, count=1)
    body = body.split("## Figures")[0] + "## Figures\n\n" + "\n".join(f"![{Path(f).stem}]({figs.name}/{f})" for f in kept) + "\n"
    import time
    out.write_text(HEAD.format(today=time.strftime("%Y-%m-%d"), run=run.relative_to(ROOT), channel=summary["channel"]) + body)
    print(f"wrote {out} ({len(body.splitlines())} lines of body, {len(kept)} figures)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
