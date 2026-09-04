# Architecture

## The pipeline as data contracts

```
ModelSpec ──realize()──► Prediction{truth: TruthTable | xsec_vector}
                                │
        ChannelSpec ────────────┼─► is_signal / in_phase_space / observables / binning
                                │
            ┌───────────────────┴──────────────────┐
   unfolded │                                       │ folded
            ▼                                       ▼
 xsec_vector_from_truth()                 expected_true_cells()  ──► Surrogate.fold() + background()
 (cm²/GeV²/nucleon per cell)              (events per true cell at the data POT)
            │                                       │
            ▼                                       ▼
 PaperRelease.compare()                    data_reco_cells() + poisson_gof()
 (MINERvA benchmark engine: total/shape χ², α, norm offset)
```

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
`save/load` with provenance.

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
parallel processes through a genie-agent environment snapshot, then `gntpc -f gst`. Target-mix
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
Save a `TruthTable` (`TruthTable.save`) or point an `external` model at a GENIE `gst` file. Quote
`sigma_per_nucleon_cm2` (the flux-averaged total for the events in the file) to enable absolute
comparisons. A new generator needs only an adapter that fills the required columns.

## Adding a channel
1. Write `channels/<name>.yaml` (copy the MINERvA one; mark every physics field's status).
2. If the observables are new, add functions to `observables.py`.
3. If the data release is not a linearised-2D MINERvA-style release, add a manifest kind to the
   benchmark bridge (`ndp/compare/minerva_bridge.py`) — see the low-recoil draft channel.
4. Build the surrogate from the experiment's paired MC (`ndp surrogate build`), check closure.

## Building the MC caches (MINERvA)
```python
from ndp.config import load_site_config; from ndp.adapters import minerva_anatuple as mad
cfg = load_site_config(); mc = cfg.data_dir / "MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root"
mad.read_truth(mc).save(cfg.data_dir / "cache/truth_mc110040.npz")
r = mad.read_reco(mc, is_mc=True)   # -> reco_mc110040.npz (+ _truthcols.npz); data -> reco_data10066.npz
```
(`ndp/cli.py::_cmd_surrogate_build` lists the exact keys it expects.)
