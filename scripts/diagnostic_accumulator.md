# scripts/diagnostic_accumulator.py

Reusable per-stage histogram accumulator for DiagnosticPayload emission.

## Purpose

Provides the producer-side companion to `scripts/diagnostic_schema.py`.
Analysis scripts import `StageAccumulator`, call `fill_*` during their
event loops, and get a validated DiagnosticPayload dict — no manual
bookkeeping of per-stage numpy arrays.

Experiment-agnostic: works for MINERvA, DUNE, T2K, or any analysis that
applies a sequential cut chain.

## API

### Construction

```python
from scripts.diagnostic_accumulator import StageAccumulator

acc = StageAccumulator(
    bin_edges=[0.0, 0.5, 1.0, 1.5, 2.0],
    stage_names=["no_cuts", "after_FV", "after_PID", "after_all_cuts"],
    cut_labels=["FV", "PID", "AllCuts"],
)
```

- `bin_edges`: N+1 values for N bins, strictly increasing.
- `stage_names`: ordered list of stage names. Stage 0 is typically
  `"no_cuts"` (all events before any selection).
- `cut_labels` (optional): N_stages - 1 labels. Stage `i >= 1` gets
  `cumulative_cuts = cut_labels[:i]`. Stage 0 gets `[]`.

### Filling (inside event loops)

| Method | Use case |
|--------|----------|
| `fill_mc(stage_index, bin_index, is_signal)` | Single MC event, single stage |
| `fill_mc_through(max_stage, bin_index, is_signal)` | MC event into stages 0..max_stage (reco-binned, reco-selected) |
| `count_mc_through(max_stage, is_signal)` | Increment scalar MC counters for stages 0..max_stage |
| `fill_eff_numerator_through(max_stage, true_bin_index)` | Efficiency numerator into stages 0..max_stage (truth-binned) — see below |
| `fill_data(stage_index, bin_index)` | Single data event, single stage |
| `fill_data_through(max_stage, bin_index)` | Data event into stages 0..max_stage |
| `count_data_through(max_stage)` | Increment scalar data counters for stages 0..max_stage |
| `fill_truth(bin_index)` | Truth event (no stage dimension) |
| `set_truth_total(n_total)` | Total truth entries processed |
| `attach_cutflow_crosscheck(external_totals, *, source, compare, rel_tol)` | Cross-check integrated stage totals vs an external C++ cutflow (v1.3+) — see below |

Out-of-range `bin_index` values (< 0 or >= n_bins) are silently ignored.

### Efficiency numerator (v1.2+)

`fill_mc_through` records the **reco-binned, reco-selected** signal — the
right histogram for purity and data/MC, but **not** an efficiency numerator.
The efficiency numerator must be matched to the `fill_truth` denominator:
events that are **truth-signal-in-phase-space AND reco-selected**, binned by
the **TRUE** kinematic variable. Fill it separately:

```python
# inside the MC reco loop, for each event:
in_ps, true_pT, _ = truth_signal_phase_space(event)   # analysis-specific
if in_ps:
    acc_pT.fill_eff_numerator_through(n_passed, find_bin(true_pT, PT_EDGES))
```

`max_stage` is the number of reco cuts the event passes; the event is recorded
in stages 0..max_stage at its **true** bin. Out-of-range true bins are dropped
(consistent with the truth denominator). Because the numerator is a subset of
the truth denominator and binned the same way, `eff_numerator / truth ≤ 1`.
Omitting it leaves the field empty and the consumer falls back to the
(unphysical) reco-binned numerator, tripping the `efficiency >= 1.0` flag.

### Cut-variable distributions (v1.1+)

Declare and fill the distribution of each variable the analysis cuts on
(observable-independent — vertex Z, apothem radius, a BDT score, a
boolean pass/fail, ...). These populate the payload's `cut_variables`
section and are rendered once by the consumer.

```python
acc.register_cut_variable(
    "vertexZ", edges=[5980.0, 6000.0, ..., 8420.0],
    label="Vertex Z", units="mm", stage="entering_ZRange")
# boolean cut: 2-bin edges
acc.register_cut_variable(
    "hasMatch", edges=[-0.5, 0.5, 1.5],
    label="Has MINOS match", units="", stage="entering_HasMatch")

# inside the event loops:
acc.fill_cut_variable("vertexZ", z_value, "data")        # kind = data
acc.fill_cut_variable("vertexZ", z_value, "mc_signal")   # or mc_signal
acc.fill_cut_variable("vertexZ", z_value, "mc_background")
```

| Method | Use case |
|--------|----------|
| `register_cut_variable(name, edges, *, label, units, stage)` | Declare a cut-variable histogram |
| `fill_cut_variable(name, value, kind)` | Fill one entry; `kind ∈ {mc_signal, mc_background, data}` |

`value` is binned against the variable's own `edges`; out-of-range
values (under/overflow) are dropped. `stage` is a free-text label
recording which selection stage the histogram was filled at (the fill
stage is the analysis's choice; for MINERvA we use "entering that cut" —
all prior cumulative cuts applied, this cut not yet — see
the analysis decision log).

### Cutflow cross-check vs an external C++ table (v1.3+)

When the selection is driven by a compiled (C++) event loop, MAT's
`runEventLoop` dumps its `mycuts` cut-summary table to
`cutflow_summary_*.json`. The Python companion pass — which re-implements the
same selection to build the per-bin histograms — cross-checks its own scalar
stage totals (from `count_mc_through` / `count_data_through`) against that
external table:

```python
import json
external = json.load(open(run_dir / "cutflow_summary_pT.json"))["stages"]
# external is a list of {"name": <stage>, "n_data": ..., "n_mc": ...}
acc.attach_cutflow_crosscheck(
    external, source="cutflow_summary_pT.json",
    compare=("n_data",),   # data counts are unweighted -> must match exactly
    rel_tol=1e-6)
```

`external_totals` may be that list (each entry carrying `name`/`stage`) or a
dict keyed by stage name. Only the quantities named in `compare` are diffed;
`"n_data"` is the safe default (unweighted integer counts). Add `"n_mc"` only
when the companion pass reproduces the C++ CV weighting, else a
weighted-vs-count difference would false-positive. Stage names must match the
accumulator's `stage_names`; a compared stage present in one source but not the
other is itself recorded as a mismatch. A resulting `status: "fail"` is
escalated by `diagnostic_plots.py` to a hard sanity violation — it means the two
selection implementations diverge, a correctness bug.

### Emission

```python
payload_dict = acc.to_payload(
    experiment="MINERvA", channel="CC-inclusive",
    variable="pT", variable_label="Muon p_T", units="GeV/c",
    signal_definition="numu CC (mc_incoming==14, mc_current==1)",
    pot_data=2.47e18, pot_mc=4.97e19,
)
```

Returns a dict conforming to `DiagnosticPayload` schema v1.3, ready for
`json.dumps` (the `cut_variables` section is included, empty if none
were registered; each stage carries `eff_numerator`; a `cutflow_crosscheck`
section is included only if `attach_cutflow_crosscheck` was called). Validate
the written file with `scripts.diagnostic_schema.load_payload()`.

## Typical usage pattern

```python
# 1. Create one accumulator per kinematic variable
acc_pT = StageAccumulator(bin_edges=PT_EDGES, stage_names=STAGES, ...)
acc_pll = StageAccumulator(bin_edges=PL_EDGES, stage_names=STAGES, ...)

# 2. MC reco loop: determine how many cuts each event passes
for event in mc_reco:
    n_passed = cuts_passed_count(event)
    bin_pT = find_bin(event.pT, PT_EDGES)          # reco bin
    acc_pT.fill_mc_through(n_passed, bin_pT, is_signal)
    acc_pT.count_mc_through(n_passed, is_signal)
    # efficiency numerator: truth-binned, truth-phase-space-gated
    in_ps, true_pT, _ = truth_signal_phase_space(event)
    if in_ps:
        acc_pT.fill_eff_numerator_through(n_passed, find_bin(true_pT, PT_EDGES))

# 3. Data loop: same structure
# 4. Truth loop: acc_pT.fill_truth(...)
# 5. Emit: acc_pT.to_payload(...)
```

## Dependencies

- numpy

## Consumers

- The vectorized per-stage CC-inclusive companion pass bundled with the
  `/ndp-execute` skill — first user (MINERvA CC-inclusive)
- Any future analysis script following the per-stage convention in
  `docs/cc_inclusive_reco_template.md`
