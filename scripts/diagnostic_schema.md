# scripts/diagnostic_schema.py

Standardized JSON schema for analysis cross-check diagnostic payloads.

## Purpose

Defines the data contract between analysis scripts (producers) and
diagnostic tools (consumers). A single validated payload file carries
everything needed to produce efficiency, purity, data/MC comparison,
and cutflow diagnostics — without the consumer knowing anything about
ROOT files, branch names, or detector geometry.

## Schema version

Current: `1.3`. Supported on load: `1.0`, `1.1`, `1.2`, `1.3`.

`1.1` adds the optional top-level `cut_variables` section (cut-variable
distributions, data vs MC). `1.2` adds the optional per-stage
`eff_numerator` field (truth-binned efficiency numerator). `1.3` adds the
optional top-level `cutflow_crosscheck` section (Python-vs-external
integrated stage-total comparison — see below). All are backward
compatible: `1.0`/`1.1`/`1.2` payloads load unchanged, the newer sections
default to `[]` / `None`.

## Structure

```
DiagnosticPayload
├── schema_version: str ("1.0")
├── metadata
│   ├── experiment: str (e.g. "MINERvA", "DUNE", "MicroBooNE")
│   ├── channel: str (e.g. "CC-inclusive")
│   ├── variable: str (e.g. "pT")
│   ├── variable_label: str (e.g. "Muon p_T")
│   ├── units: str (e.g. "GeV/c")
│   ├── signal_definition: str (free text)
│   ├── pot_data: float (> 0)
│   └── pot_mc: float (> 0)
├── binning
│   └── edges: list[float] (strictly increasing, N+1 values for N bins)
├── truth
│   ├── counts: list[float] (N bins, generator-level signal)
│   └── n_total: int (total truth entries processed)
├── stages: list[Stage] (≥ 1)
│   └── Stage
│       ├── name: str (e.g. "after_fiducial", "after_all_cuts")
│       ├── cumulative_cuts: list[str] (cut labels applied so far)
│       ├── mc_signal: list[float] (N bins, raw MC counts)
│       ├── mc_background: list[float] (N bins, raw MC counts)
│       ├── data: list[float] (N bins, raw data counts)
│       ├── n_mc_signal: int
│       ├── n_mc_background: int
│       ├── n_data: int
│       └── eff_numerator: list[float] | None (optional, v1.2+; N bins,
│              truth-binned efficiency numerator — see below)
├── cut_variables: list[CutVariable] (optional, default []; v1.1+)
│   └── CutVariable
│       ├── name: str (e.g. "vertexZ"; used in the dist_{name}.pdf filename)
│       ├── label: str (axis label, e.g. "Vertex Z")
│       ├── units: str (e.g. "mm"; "" for dimensionless / boolean)
│       ├── stage: str (selection stage histogram was filled at)
│       ├── edges: list[float] (M+1 strictly increasing; use 2-bin edges
│       │                       like [-0.5, 0.5, 1.5] for boolean cuts)
│       ├── mc_signal: list[float] (M bins, raw MC counts)
│       ├── mc_background: list[float] (M bins, raw MC counts)
│       └── data: list[float] (M bins, raw data counts)
└── cutflow_crosscheck: CutflowCrosscheck | None (optional, v1.3+)
    ├── source: str (provenance, e.g. the cutflow_summary_*.json path)
    ├── status: str ("pass" | "fail"; must match the comparisons)
    ├── rel_tol: float | None (match tolerance used)
    └── comparisons: list[dict] (each: stage, quantity, python, external,
                     rel_diff, match; a false `match` => status "fail")
```

### Efficiency numerator (`eff_numerator`, v1.2+)

`mc_signal` is the **reco-binned, reco-selected** signal (the framework's
`selectedSignalReco`) — correct for purity and data/MC, but **not** a valid
efficiency numerator. The efficiency numerator must be the population matched
to the `truth.counts` denominator: MC events that are **truth-signal-in-phase-
space AND pass the reco selection**, binned by the **TRUE** kinematic variable
(the framework's `efficiencyNumerator`). `eff_numerator` carries exactly that,
per stage. Because it is a subset of the truth denominator and binned the same
way, `eff_numerator / truth.counts ≤ 1` — a probability, as efficiency must be.

When absent (1.0/1.1 payloads), the consumer falls back to `mc_signal` for the
efficiency (the original behaviour), and the `efficiency >= 1.0` sanity flag in
`diagnostic_plots` fires, signalling that the payload predates the fix. New
producers should always emit `eff_numerator`.

`cut_variables` carry the distribution of each variable the analysis
cuts on. They are **observable-independent** (the same vertex-Z or BDT
distribution applies whether you bin the result in pT or p∥), so an
analysis emits the same `cut_variables` block in every per-observable
payload. The consumer renders them once (see `diagnostic_plots.md`).
Each cut variable has its own binning, independent of `binning.edges`.

### Cutflow cross-check (`cutflow_crosscheck`, v1.3+)

When a compiled (C++) event loop drives the selection, MAT's `runEventLoop`
prints a `mycuts` cut-summary table (per-stage counts / % Eff / % Purity for
MC, total entries for data) and dumps it to `cutflow_summary_*.json`. The
Python companion diagnostic pass — which independently re-implements the same
selection to build the per-bin, per-stage histograms — cross-checks its own
integrated stage totals against that external table via
`StageAccumulator.attach_cutflow_crosscheck(...)`. A `status: "fail"` means the
two selection implementations diverge (a correctness bug, not a fluctuation);
`diagnostic_plots.py` escalates it to a hard sanity violation so the run cannot
be quoted. Absent (`None`) for pure-Python loops that need no such cross-check.

## Validation

`load_payload(path)` reads JSON, parses into dataclasses, and runs
`validate_payload()` which checks:

- Schema version is one of the supported versions (`1.0`, `1.1`, `1.2`, `1.3`)
- Bin edges are strictly increasing with ≥ 2 edges
- All per-bin arrays have length matching `binning.n_bins`
- Each stage's `mc_signal`/`mc_background`/`data` bins are finite,
  non-negative *numbers* (rejects NaN/inf, negatives, and non-numeric
  dtypes including `bool`/`str`)
- All event counts are ≥ 0
- `eff_numerator` (if present) has length `binning.n_bins` and is
  non-negative
- POT values are > 0
- At least one stage is present
- Stage `name`s are unique, and a `"no_cuts"` stage is present (the
  pre-cut baseline). The two legacy v1.0 single-stage `after_all_cuts`
  payloads predate this convention and are out of the certified corpus.
- Each `cut_variables` entry (if present) has strictly increasing edges
  and `mc_signal`/`mc_background`/`data` lengths matching its bin count,
  with non-negative counts
- `cutflow_crosscheck` (if present) has a `status` of `"pass"`/`"fail"`,
  each comparison carries `stage`/`quantity`/`match`, and the `status` is
  consistent with the comparisons (fail iff any `match` is false)

Raises `PayloadError` (subclass of `ValueError`) on any violation.

## Usage

```python
from scripts.diagnostic_schema import load_payload

payload = load_payload("path/to/payload.json")
# payload.metadata.experiment → "MINERvA"
# payload.binning.n_bins → 14
# payload.stages[0].mc_signal → [12.0, 45.0, ...]
```

## Conventions

- All MC counts are **raw** (un-scaled). The consumer applies
  `pot_data / pot_mc` scaling.
- `truth.counts` is the generator-level signal distribution in the same
  binning. Used as the efficiency denominator.
- `stages` are ordered by increasing restrictiveness (more cuts applied).
  The last stage is typically the final selection.
- `cumulative_cuts` lists all cut labels applied up to and including
  this stage, in application order.

## Best practices

- **Always emit all cut stages**, not just the final selection. This
  enables cut-evolution diagnostics (efficiency sculpting, purity
  improvement, background rejection through the chain).
- **Stage naming:** `"no_cuts"` for the pre-selection baseline, then
  `"after_{CutLabel}"` for each cumulative cut (e.g. `"after_ZRange"`,
  `"after_Apothem"`, ..., `"after_all_cuts"` or `"after_{LastCut}"`).
- **`cumulative_cuts`** should list all cut labels applied so far, in
  application order. Stage 0 gets `[]`; stage `i` gets `labels[:i]`.
- Use `scripts/diagnostic_accumulator.StageAccumulator` as the standard
  producer. It handles the per-stage numpy bookkeeping and emits a
  validated payload dict.

## Consumers

- `scripts/diagnostic_plots.py` — standardized diagnostic plots + summary JSON

## Producers

- `scripts/diagnostic_accumulator.py` — reusable per-stage histogram accumulator
