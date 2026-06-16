# Systematic uncertainties — status summary (2026-06-16)

Consolidated status of the systematics stage (M4), mapped onto the paper's
categories. The paper validates against **4 ancillary covariance files**
(`cov_stat`, `cov_flux`, `cov_energyscale`, `cov_total`); everything else is
checkable only via Fig 8 (the per-category fractional-uncertainty breakdown) or
the assembled total.

Paper §9 names **three categories** — Flux, Detector, Interaction model — while
**Fig 8** plots seven curves that refine them. Status by Fig 8 curve:

| Fig 8 curve | paper category | our status | ours (1A, median/cell) | target |
|---|---|---|---|---|
| **Muon Reconstruction** | Detector | ✅ done | **3.79 %** | dominant; Fig 8 bulk 3–5 % |
| Flux | Flux | ✅ normalization; shape deferred | 3.23 % | anc `cov_flux` ~4.0 % |
| Statistical | — | ✅ toys (1A) | 4.70 % | anc `cov_stat` (scaled) 0.90 |
| Models (interaction) | Interaction model | ◑ GENIE+2p2h+RPA; GenieRvx1pi missing | 1.2 % | Fig 8 ~2–3 % |
| Normalization | Flux | ❌ not built | — | ~1.4 % |
| Hadronic Response | Detector | ❌ not built | — | ~1 % |

## Validated against dedicated anc files

| group | ours | anc | metric |
|---|---|---|---|
| Flux normalization | 3.23 % | 4.01 % | within ~19 % (norm model; shape deferred) |
| Statistical (1A toys) | 4.70 % | 5.33 % (scaled to 1A) | 0.90 |
| **Muon energy scale** (MINERvA ⊕ MINOS) | 3.07 % | 3.55 % | **median per-cell ratio 0.99** (was 0.84) |

### Energy-scale caveat (median ratio ≠ magnitude ratio)
The 0.99 is the **median of the per-cell ratios** — the typical cell matches.
But the **magnitude** ratio is median(ours)/median(anc) = 0.87, and the residual
is **structured**: we slightly over-predict the small-uncertainty bulk cells and
**under-predict the high-uncertainty edge cells** (steep muon-momentum-peak
gradient) by ~18 %. Likely cause: the deferred MINOS-band flux-weight
correlation, largest at the spectrum edges. See
`img/muon_1A/energyscale_distributions.png`.

## Folded into Cov_total only (no dedicated anc file)

| group | ours (1A, median/cell) | note |
|---|---|---|
| GENIE (56 universes) | 1.19 % | small for inclusive (knobs cancel in ratios) |
| 2p2h low-recoil tune (3 universes) | 0.31 % | near-pure shape; integrated σ unchanged |
| RPA (4 universes: HighQ2/LowQ2) | 0.10 % | negligible for inclusive; integrated σ unchanged |
| MINOS efficiency | 1.37 % | exact `GetWeightRatioToCV` |
| Beam angle X / Y | 0.34 / 0.30 % | re-applies the 20° cut |
| Muon resolution / angle resolution | ~0.1 % | negligible, as predicted |

## Cov_total (playlist 1A)

| | systematic-only | total |
|---|---|---|
| ours | **5.30 %** | 7.30 % |
| published anc | 6.26 % | 6.83 % |
| ratio | **0.85** | 1.05 (overshoots) |

Systematic-only is **85 %** of the anc systematic budget. The total overshoots
only because the 1A statistical band (4.7 %) stands in for the much smaller
full-dataset stat (~1.5 %) — the **12-playlist combine** resolves this and is
the gate for a clean `cov_total` validation.

## Covariance conventions (from the MAT source)

- Many-universe band (flux 1000, 2p2h 3): **sample covariance about the universe
  mean, 1/N** (`MnvVertErrorBand::CalcCovMx`, `fUseSpreadError=false` for >1
  universe).
- ±1σ pair (GENIE knob, each muon band): `outer((σ⁺−σ⁻)/2)` — the N=2 special
  case of the above (verified).
- Total = Σ group covariances + Cov_stat.

## Remaining

1. **GenieRvx1pi** — non-resonant single-π normalization (confirmed a separate
   band in `GetStandardSystematics`); the likely Fig 8 "Models" residual.
2. **geant4 + response** — Geant-hadron + calorimetric response (the Fig 8
   "Hadronic Response" curve, ~1 %), and a **Normalization** band (~1.4 %).
3. **Flux shape** term (per-cell de-correlation; PPFX universes already loaded).
4. **MINOS-band flux-weight correlation** (closes the energy-scale edge deficit).
5. **Full-12-playlist combine** — shrinks 1A stat 4.7 % → ~1.4 %; the gate for
   validating `cov_total` cleanly.

## Per-stage detail
`2026-06-16_{flux,stat,energyscale,total,2p2h,muon_reconstruction}_systematic.md`.
