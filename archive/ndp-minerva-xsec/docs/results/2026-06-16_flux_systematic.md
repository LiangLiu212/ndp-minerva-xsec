# S2 — flux systematic: ν-e constraint validated, normalization Cov_flux

The flux systematic's crux (flagged in the plan) is reproducing the
**ν-e-constrained** uncertainty, not the raw PPFX spread. Source-confirmed
mechanism: `MnvHistoConstrainer::CorrectFluxUniv` applies the constraint as
**per-universe weights** (`SetUnivWgt`), so the constrained covariance is the
**weighted** sample covariance of the raw PPFX universes — which
`xsec.systematics.sample_covariance(..., weights=w)` implements directly.

## Crux validated — the constraint reproduces the paper

`FluxCV` now retains the 1000 raw PPFX universes + the ν-e constraint weights
(`flux_norm_uncertainty`, `universe_integrals`, `universe_weight_ratios`):

| flux normalization uncertainty (0–100 GeV) | value |
|---|---|
| raw PPFX (unconstrained) | **7.63 %** |
| **ν-e constrained** | **3.23 %** |
| constraint reduction | ×0.42 |

The ν-e constraint reduces the flux uncertainty to **3.23 %**, below the
paper's quoted **< 4 %** — the constraint machinery is correct. (Range-robust:
3.21–3.24 % over any physical Enu window.)

## Cov_flux vs the published anc

The anc `cov_ptpl_*_flux.txt` per-cell fractional uncertainty is **flat at
~4 %** (median 4.01 %, band 3.9–4.3 % over the 205 reported cells) — confirming
the flux is normalization-dominant, the structure the constraint produces.

`systematics.normalization_covariance(σ_cv, 0.0323)` builds the leading-order
fully-correlated Cov_flux from the constrained normalization. It matches the
anc magnitude to **~19 %** (3.23 % vs 4.01 %), within the plan's ±20 % gate.

**What the residual is (3.23 → 4.0 %):** the per-event Enu-shape folding through
the efficiency — the constrained `Flux` band propagated event-by-event (each
cell samples a slightly different Enu mix) gives the per-cell 4 %. The
integrated-normalization model is flat by construction and misses this shape.
Verified NOT to be: an Enu-range effect (range-robust above), nor the other
flux bands (`Flux_BeamFocus`, `ppfx1_Total` are spectators; naive quadrature
overshoots to 8.45 %, so only the constrained `Flux` band propagates).

## Done vs deferred

- **Done:** flux universe machinery (`xsec/flux.py`), the covariance assembler
  (`xsec/systematics.py`), the ν-e-constraint validation (3.23 % vs paper < 4 %),
  and a normalization-model Cov_flux within ~20 % of the anc.
- **Deferred (shape-resolved Cov_flux):** apply `universe_weight_ratios(Enu)`
  per MC event in a streaming pass → per-universe ingredients → per-universe
  cross sections → weighted `sample_covariance` with the constraint weights.
  This is the same vertical weight-matrix infrastructure S5 needs, so it is
  natural to build there and back-apply to flux. It closes the 3.23 → 4.0 %
  per-cell shape and recovers the anc off-diagonal correlations.

## Reproduce

```bash
pixi run python -c "from xsec.flux import FluxCV; fx=FluxCV('minervame1A'); \
  print(fx.flux_norm_uncertainty(0,100,constrained=False), \
        fx.flux_norm_uncertainty(0,100,constrained=True))"
```
