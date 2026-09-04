# 2p2h low-recoil tune systematic (interaction model)

The MnvTune 2p2h band has **three universes** (`Get2p2hSystematics`,
MnvTuneSystematics.cxx:64–73), each at nsigma=1 with a different MEC-pair fit:

| universe | variation | acts on |
|---|---|---|
| 1 | nn/pp-only fit | MEC, struck nn or pp pair (`mc_targetNucleon-2000000200` ∈ {0,2}) |
| 2 | np-only fit | MEC, struck np pair (== 1) |
| 3 | QE→2p2h fit | true CCQE (`mc_intType==1`) on nuclei |

2p2h is **vertical** (weight-only): the reco/true slots are the CV ones; each
universe swaps the CV low-recoil weight for its variation
(`weights.twop2h_variation_ratio` = variation / CV, applied on top of the full
CV stack). `make_2p2h_universes.py` streams the MC reco + Truth trees once
(RunLog 2026_06_16_182948, 41 files, **0 failures, 468 s, 8 workers**) and fills
the three per-universe ingredient sets (migration, bkg, eff_denom); eff_num is
each migration's column sum.

## Covariance convention (read from the MAT source)

The decisive point is how MINERvA turns 3 universes into a covariance. A vertical
band with >1 universe has `fUseSpreadError = false`
(MnvVertErrorBand.cxx:52–55), so `CalcCovMx` (lines 358–392) is the **sample
covariance about the universe *mean*, normalized by 1/N**:

  Cov[i,k] = (1/3) Σ_u (σ_u[i] − σ̄[i])(σ_u[k] − σ̄[k]),  σ̄ = (σ_1+σ_2+σ_3)/3

i.e. `xsec.systematics.sample_covariance` — **not** a ±1σ pair. (Cross-check:
the 2-universe GENIE `pair_covariance` is exactly this formula for N=2, so the
framework is internally consistent.) Consequence worth stating: the band
captures the **spread among the three model variations**, not their offset from
the CV tune.

## Result (playlist 1A)

| metric | value |
|---|---|
| 2p2h fractional uncertainty / cell (median, reported cells) | **0.31 %** |
| p16–p84 band | [0.03, 1.09] % |
| max (sparse high-p_T / low-p_∥ corner) | 8.1 % |
| per-universe median \|σ_u − σ̄\|/σ (nn/pp, np, QE) | 0.24 % / 0.18 % / 0.22 % |
| **integrated** σ across the 3 universes | 1.0001 / 1.0001 / 0.9998 |
| cov_2p2h symmetric / PSD | yes / yes (rank 2, min eig ≈ −1e-95) |

**2p2h is almost purely a shape effect** — the integrated inclusive rate is
data-driven and barely moves (≤0.03 %), while the band redistributes within the
(p_T, p_∥) plane at the few-tenths-of-a-percent level, reaching ~8 % only in a
statistics-starved corner cell. This is the expected size for an inclusive
muon-kinematics measurement and matches the paper's statement that the
interaction model is "not dominant in any (p_T, p_∥) bin."

It is small because 2p2h enters the cross section only through MC ratios that
largely cancel: the migration is conditional (row-renormalized in unfolding),
and the efficiency num/denom both carry the same reweight.

## Effect on Cov_total (1A)

| | systematic-only (median/cell) | total (median/cell) |
|---|---|---|
| flux ⊕ es ⊕ GENIE | 4.70 % | 6.93 % |
| **+ 2p2h** | **4.78 %** | **6.93 %** |
| published anc | 6.26 % (systematic) | 6.83 % (total) |

2p2h lifts the systematic-only budget 4.70 → 4.78 % (now **76 %** of the anc
systematic) and is negligible in the total, which is still dominated by the 1A
statistical band (4.7 %). No dedicated anc file exists for 2p2h (it folds into
`cov_total` with the rest of the interaction model).

## What this leaves (interaction model)

GENIE knobs (done, 1.19 %) and 2p2h (done, 0.31 %) are in. Remaining in this
category: **RPA** Valencia high-Q²/low-Q² variations (the CV `RPAWeight` is
already ported to parity 1e-12; only the variation modes need wiring) and a
check of whether the non-resonant single-π reduction needs its own band or is
covered by the GENIE Rvn/Rvp knobs.

## Reproduce

```bash
pixi run python make_2p2h_universes.py --workers 8 --playlist minervame1A
pixi run python assemble_total.py --ingredients <ing> --xsec <xsec> \
    --energyscale <es> --stat-cov <stat> --genie <genie> --twop2h <twop2h>
```
