# Architecture

## The workflow: forward folding

Three ingredients, three roles:

| ingredient | role | where |
|---|---|---|
| **experiment data** (MINERvA open-data AnaTuples; the published cross section where one exists) | what reality recorded: selected reconstructed candidates in *any* reco observable the analyst defines; optionally the experiment's unfolded result on its own grid | `ndp/adapters/`, `ndp/compare/` |
| **experiment official MC** (paired truth ↔ reco) | the detector's performance: what the reconstruction and selection do to true signal. It is used only to *learn a surrogate*; it is never the "model" being tested (except as the `reference_mc` control) | `ndp/surrogate/` |
| **theorist model** (generator run, reweight, external sample) | the physics hypothesis, as truth-level events | `ndp/theory/` |

The model's truth is pushed **forward** through the surrogate into reconstructed space and compared
with the data counts there. Nothing is unfolded by the platform: the experiment's unfolded release is
used only as a second, independent comparison where the measurement is the one they published.

```
ModelSpec ──realize()──► Prediction{truth: TruthTable | xsec_vector}
                                │
        ChannelSpec ────────────┼─► is_signal / in_phase_space / selection / data files / normalisation
        Measurement ────────────┼─► (x, y) truth observables ↔ (x, y) reco observables, edges, [release]
                                │
            ┌───────────────────┴──────────────────┐
     folded │  (primary; any measurement)           │ unfolded  (only where the experiment published this grid)
            ▼                                       ▼
 expected_true_cells()  ──► Surrogate.fold()    xsec_vector_from_truth()
 (events per true cell at the data POT)         (cm²/GeV²/nucleon per cell)
            │  + background()                       │
            ▼                                       ▼
 data_reco_cells() + poisson_gof()              PaperRelease.compare()
 (−2lnL, Pearson with MC stat, data/pred)       (MINERvA benchmark engine: total/shape χ², α, norm offset)
```

### Measurement (`measurements/<channel>/*.yaml`, `ndp/channels/measurements.py`)
The analyst's choice of observables. Each axis has a **truth** side (a name from
`ndp/channels/observables.py` or an expression over truth columns + observables, e.g. `"E_nu - lep_E"`)
and a **reco** side (a name from `ndp/channels/reco_observables.py` or an expression over the cached
reco columns, e.g. `"reco_E_mu + reco_recoil_E"`), plus edges, units and a plotting flag. Omit `y` for a
one-dimensional measurement. Every channel has a `published` measurement (its `binning`, with the
paper's release attached); the rest are folded-only. A measurement is compared with data through a
surrogate learned *for that measurement* (`ndp surrogate build --measurement <name>`), stored under
`surrogates/<channel>/<measurement>/`, and certified by closure on the training MC.

### Reco cache (`ndp data cache`)
The adapter caches, per AnaTuple, the selection result and the reconstructed quantities user
observables are built from (`RECO_CACHE_COLUMNS` in `adapters/minerva_anatuple.py`: muon p/θ/pT/p∥/E,
MINOS momentum, the recoil-energy family, the analysis tool's E_ν/W/x/y, visible energy, track
multiplicities, vertex). Caches carry a version; a stale cache is refused with the rebuild command.

### TruthTable (`ndp/events.py`)
Columnar truth events in GeV / GeV² / mm: neutrino, primary lepton 4-momentum, current,
interaction type (NDP codes 1 QE, 2 RES, 3 DIS, 4 COH, 5 MEC, 0 other), target Z/A, Q², W,
weight, optional vertex, optional jagged final-state particles (`fs_offsets` + `fs_*`).
`meta["norm"]` says what a unit of weight means:

| kind | meaning | who produces it |
|---|---|---|
| `xsec_per_nucleon` | σ_cell = Σw · `xsec_per_unit_weight` (cm²/nucleon) | GENIE runs (flux-averaged σ from splines), external samples with a quoted σ |
| `pot` | the sample represents `pot` protons on target on the modelled target | the experiment's own MC |
| `shape` | no absolute scale | external samples without a quoted σ (unusable for rates) |

`meta["has_geometry"]` tells the channel whether the fiducial-vertex phase-space cut applies
(the experiment's truth: yes; a point-target generator run: no — its target mass is `n_nucleons`).

### ChannelSpec (`channels/*.yaml`, `ndp/channels/registry.py`)
Everything about the measurement that is not the theorist's model. Physics fields carry a
`status`: `decided` (ratified by the physicist / taken from the publication), `default` (platform
choice with recorded evidence, awaiting ratification), `open` (must be settled before use).
Observables are named functions in `ndp/channels/observables.py`; the binning is a linearised 2D
grid whose cell formula is written exactly as the paper writes it (`ipt*n_pz + ipz`).

### Surrogates (`ndp/surrogate/`)
`Surrogate.fold(true_cells) -> reco_cells`, `fold_events(x, y, w)`, `background(pot)`, plus
`save/load` with provenance. One surrogate per (channel, measurement): `ndp/surrogate/build.py`
assembles the training arrays from the cached MC (denominator from the Truth tree, numerator +
migration from the reco rows' truth branches, background from every other selected candidate) for
the measurement's observables and runs the closure check.

- **BinnedResponse**: `reco = P @ (eff · true) + bkg·POT`. `eff[j] = num[j]/den[j]` (den = signal in
  the true phase space from the Truth tree; num = reco-selected signal in the phase space), `P[:, j] =
  M[:, j]/num[j]` so columns sum to ≤ 1 (loss to outside the reco grid is a loss), and `bkg` = every
  other selected event in the reco grid (non-signal + signal whose truth is outside the phase space)
  per POT. Folding the training truth reproduces the training reco cell by cell — the closure test.
- **SmearingSurrogate**: per true cell and axis, robust location/width of reco−true (pT) and
  reco/true (p∥), plus the acceptance; samples reconstructed events from any truth sample, so it
  works on any binning. Gaussian core only — tails are not modelled (documented limitation). Its
  interface is what a conditional normalising flow or diffusion model would implement.

### Normalisation constants (channel `normalization`)
`phi_per_pot_cm2` (integrated flux, 0–100 GeV) and `n_nucleons` (fiducial target) convert a
POT-normalised MC into a per-nucleon cross section and an absolutely normalised model into an
event rate. Both are the published values the MINERvA 2D reproduction adopted. The published
2110.13372 flux table (`resources/flux/`) drives GENIE and cross-checks Φ (its 0–100 GeV integral
is 6.27e-8 vs the 6.32e-8 quoted, the difference being the table's two-decimal truncation).

### GENIE runs (`ndp/theory/generator.py`)
`gevgen -f flux.root,flux -t <mass-fraction mix> --cross-sections <splines> --tune …` in `n_jobs`
parallel processes through an environment snapshot (`ndp.yaml: genie_env_json`), then
`gntpc -f gst`. The default snapshot is the in-repository GENIE built by `pixi run build-genie`
(`external/genie/Generator`, R-3_06_02, Pythia6 via ROOTEGPythia6, LHAPDF 6, conda-forge ROOT);
a snapshot of any other installation (e.g. a genie-agent `config/env/*.json`) is a drop-in. Target-mix
weights are mass fractions (GENIE's `GMCJDriver` treats them as density-weighted path lengths and
divides by A). Normalisation: σ_avg/nucleon = ∫Φ(E) Σ_i (w_i/A_i) σ_i(E) dE / ∫Φ dE with σ_i the
summed CC splines of each nuclide. Productions are cached by spec fingerprint.

## Reading the two comparison modes together
Unfolded space asks "does the model's cross section match what the experiment published?", with
the experiment's full covariance but also its unfolding model dependence. Folded space asks "would
the detector have recorded these counts?", with no unfolding but with the surrogate's (MC-derived)
response and only the statistical error of the small open-data slice. They answer different
questions; the platform reports both and never merges them into one verdict.

## Bringing your own events
Save a `TruthTable` (`TruthTable.save`) or point an `external` model at a generator file:

| `format` | reader | normalisation it derives |
|---|---|---|
| `genie_gst` | `adapters/genie_gst.py` | none (quote `sigma_per_nucleon_cm2`) |
| `nuwro_root` | `adapters/nuwro_root.py` (PyROOT + `event1.so`) | `event.weight` = σ_tot [cm²], taken per nucleon |
| `gibuu_finalevents` | `adapters/gibuu_finalevents.py` | perweight in 1e-38 cm²/nucleon, ÷ n_runs |
| `nuhepmc` | `adapters/nuhepmc.py` (pyhepmc) | G.C.4 flux-averaged total, unit + PerAtom/A → cm²/nucleon |

Process labels are mapped onto the NDP codes (QE, RES, DIS, COH, MEC); the generator's own code
is kept in an extra column (`gibuu_production_id`, `process_id`, `nuwro_dyn`). Quote
`sigma_per_nucleon_cm2` in the model spec to override a derived normalisation. A new generator
needs only an adapter that fills the required columns.

## Adding a measurement (your own observable)
1. Write `measurements/<channel>/<name>.yaml`: `x` (and optionally `y`) with `observable` (truth) and
   `reco`, edges, units. Expressions are allowed on both sides; `python -m ndp measurements --channel
   <channel>` lists what resolves.
2. If a reco quantity you need is not cached, add its branch to `RECO_EXTRA_BRANCHES` in the adapter,
   bump `RECO_CACHE_VERSION`, and rebuild with `python -m ndp data cache --channel <channel>`.
3. `python -m ndp surrogate build --channel <channel> --measurement <name>` (binned + parametric;
   closure is checked and printed).
4. `python -m ndp run models/<model>.yaml --channel <channel> --measurement <name>`.

## Adding a channel
1. Write `channels/<name>.yaml` (copy the MINERvA one; mark every physics field's status).
2. If the observables are new, add functions to `observables.py` / `reco_observables.py`.
3. If the data release is not a linearised-2D MINERvA-style release, add a manifest kind to the
   benchmark bridge (`ndp/compare/minerva_bridge.py`) — see the low-recoil draft channel.
4. Build the caches (`ndp data cache`) and the surrogate (`ndp surrogate build`), check closure.

## Building the MC caches (MINERvA)
```bash
python -m ndp data cache --channel minerva_me_cc_inclusive_ptpz          # data + MC (Truth tree ~15 s, reco ~1 min)
python -m ndp data status                                                 # cache versions
```
Files: `truth_<tag>.npz` (TruthTable of the Truth tree), `reco_<tag>.npz` (selection + reco columns,
versioned), `reco_<tag>_truthcols.npz` (truth branches of the reco rows), `<tag>` = `mc110040` / `data10066`.
