# scripts/diagnostic_plots.py

Standardized diagnostic plots from a `DiagnosticPayload` JSON file.

## Purpose

Consumes the experiment-agnostic payload defined by
`scripts/diagnostic_schema.py` and produces four core plot types plus a
machine-readable summary. Works for any neutrino cross-section
experiment (MINERvA, DUNE, T2K, ...) — no ROOT dependency.

## Plots produced

Outputs split into **per-observable** artifacts (written to
`--output-dir`) and **observable-independent** artifacts (written to
`--shared-dir`, once — see "Shared artifacts & dedup" below).

### Per-observable (one set per pT, p∥, ...)

For each **stage** in the payload:

1. **Efficiency vs variable** — `efficiency_{stage}_{variable}.pdf`
   Per-bin efficiency = `eff_numerator / truth.counts`, with binomial error
   bars. The numerator is the **truth-binned** count (truth signal in phase
   space AND reco-selected), carried by schema v1.2+.

   Legacy 1.0/1.1 payloads have no `eff_numerator`, so the reco-binned
   `mc_signal` is used instead — the original behavior, which can give
   efficiency ≥ 1 because numerator and denominator are then different
   populations. The `eff >= 1` sanity flag still fires on those, surfacing that
   the payload predates the fix rather than reporting a wrong-but-plausible
   number.

2. **Signal purity vs variable** — `purity_{stage}_{variable}.pdf`
   Per-bin purity = `mc_signal / (mc_signal + mc_background)`.

3. **Data/MC ratio** — `data_mc_{stage}_{variable}.pdf`
   Upper panel: POT-scaled MC (signal+bkg stacked) with data overlay. The
   first MC series is the **signal+background total** drawn in the
   `mc_signal` color and labeled "MC total (signal+bkg, POT-scaled)" — in
   `step` mode that curve traces the total, so the signal is the gap above
   the background curve. Lower panel: data / MC total ratio with ±10%
   shaded band.

### Observable-independent (written once, shared across observables)

4. **Cutflow** — `cutflow.pdf`
   Side-by-side bars showing data and MC event counts at each stage.
   **MC is POT-scaled to data** (`pot_data / pot_mc` from payload
   metadata — generic, no experiment-specific constant), so the data
   and MC bars are directly comparable. (Integrated counts per stage are
   identical for every observable, so the filename carries no variable
   suffix.)

5. **Cut-variable distributions** — `dist_{name}.pdf` (one per
   `cut_variables` entry; only when the payload carries them).
   Data vs MC for each variable the analysis cuts on. **MC is scaled to
   data via POT** (`pot_data / pot_mc` from payload metadata — generic,
   no experiment-specific constant), so data and MC overlay directly.
   Upper panel: POT-scaled MC (signal+bkg stacked) with data points
   (Poisson error bars). The first MC series is the **signal+background
   total** drawn in the `mc_signal` color and labeled "MC total
   (signal+bkg, POT-scaled)" — in `step` (outline) mode that curve traces
   the total, so the signal contribution is the gap above the background
   curve, not the labeled curve itself. Lower panel: **Data / MC ratio
   with ±10% shaded band** — built via the shared
   `comparison_plot.data_mc_figure` builder (the ratio-panel
   convention), the same builder the per-observable data/MC plots use.

## Shared artifacts & dedup

The cutflow and cut-variable distributions are observable-independent,
so they belong in one shared location, not duplicated per observable.

- `--shared-dir` (default: same as `--output-dir`) receives `cutflow.pdf`
  and the `dist_{name}.pdf` set.
- `--emit-shared` controls (re)writing of shared artifacts:
  - `auto` (default): write only if the target file is absent. Running
    for pT then p∥ into the same `--shared-dir` emits the shared set
    exactly once (the second run finds them present and skips).
  - `always`: always (re)write.
  - `never`: skip shared artifacts entirely (per-observable plots only).

`run_diagnostics(...)["shared_artifacts"]` lists what was actually
written this call.

## Summary output

`diagnostic_summary.json` — per-stage metrics:
- `efficiency_per_bin`, `purity_per_bin`, `data_mc_ratio_per_bin`
- `integrated_efficiency`, `mean_purity`, `chi2`, `ndf`
- `sanity_flags`: all warnings, human-readable (efficiency >= 1 per-bin or
  integrated — unphysical, signals a numerator/denominator mismatch
  upstream; purity outside [0,1]; data/MC ratio outside [0.1, 10])
- `sanity_status`: `"pass"` or `"fail"`. **`fail` means the run is not a
  quotable result.** It is set by either a hard physical-bound violation or
  a declared-range violation (below). The CLI and producing event loops
  exit non-zero on `fail`; artifacts are still written (stamped `fail`) for
  debugging.
- `hard_violations`: the subset of impossibilities that force `fail` —
  efficiency or purity outside [0, 1] (universal physical bounds, always
  checked), **and** any cutflow cross-check mismatch (below). A
  suspicious-but-possible value (e.g. a large data/MC ratio) is an advisory
  `sanity_flag` only, **not** a hard violation.
- `cutflow_crosscheck` (schema v1.3+, present only when the payload carries
  it): the Python-vs-external integrated stage-total comparison echoed from
  the payload (`source`, `status`, `rel_tol`, `comparisons`). A `status` of
  `"fail"` means the Python companion selection and the compiled C++
  selection diverge — a correctness bug — so each mismatched comparison is
  pushed into `hard_violations`, forcing `sanity_status: "fail"` at the same
  severity as an unphysical efficiency.
- `expected_range_violations`: violations of the optional analysis-declared
  `expected_ranges` passed to `run_diagnostics` (also force `fail`).
  Supported keys: `efficiency_per_bin`, `integrated_efficiency`,
  `mean_purity` (every stage); `final_integrated_efficiency`,
  `final_mean_purity` (last stage only). Each maps to an inclusive
  `(lo, hi)`. An unknown key raises (typo protection). Declaring a range up
  front turns a "looks plausible but is wrong" benchmark into an immediate
  fail rather than an ignored warning.
- `shared_artifacts`: shared files written this call (may be empty if
  already present / `--emit-shared never`)
- `cut_variables`: per-cut-variable `{name, stage, n_data, n_mc,
  data_mc_ratio_per_bin}` (n_mc POT-scaled; ratio = data / POT-scaled MC
  total, `null` where MC total is zero — recorded even when the shared
  PDF was skipped by the dedup policy)

## Invocation

```bash
# Two observables sharing one cutflow + cut-var distribution set:
python scripts/diagnostic_plots.py payload_pT.json \
    --output-dir results/pT  --shared-dir results/shared
python scripts/diagnostic_plots.py payload_pll.json \
    --output-dir results/pll --shared-dir results/shared
```

As a library:

```python
from scripts.diagnostic_plots import run_diagnostics
summary = run_diagnostics("payload.json", "results/pT",
                          shared_dir="results/shared")  # emit_shared="auto"
```

## Inputs

A JSON file conforming to the `DiagnosticPayload` schema. Versions 1.0-1.3
are accepted (`diagnostic_schema.SUPPORTED_VERSIONS`); the accumulator emits
1.3. See `scripts/diagnostic_schema.md` for the reference.

## Conventions

- MC counts in the payload are **raw**; this tool applies
  `pot_data / pot_mc` scaling internally.
- Efficiency denominator is `truth.counts` (same for all stages).
- Bins with zero truth or zero MC total produce `NaN` (shown as gaps).
- Plot style: every color/marker comes from `plot_style`; PDF output. (The
  bundled default palette is `solid-dark`, not Okabe-Ito — see
  `scripts/plot_style.md`.)

## Dependencies

- Python 3.10+
- numpy
- matplotlib
