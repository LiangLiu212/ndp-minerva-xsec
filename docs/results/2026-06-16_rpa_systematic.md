# RPA systematic (interaction model)

The Valencia RPA tune uncertainty is **two ±1σ bands → 4 universes**
(`GetRPASystematicsMap`, MnvTuneSystematics.cxx:119-127):

| band | variations | weightRPA component |
|---|---|---|
| **HighQ2** | 1 (+), 2 (−) | enhancement: CV ± 0.6·(non-rel − CV), with traps |
| **LowQ2** | 3 (+), 4 (−) | suppression: CV ± 0.25·(1 − CV) when CV<1 (muon-capture) |

RPA is **vertical** (weight-only), gated on true CCQE (`mc_intType==1`) on Z≥6;
each universe swaps the CV RPA weight for its band variation
(`weights.RPAWeight.variation_ratio`). `make_rpa_universes.py` runs the 4-universe
vertical pass (RunLog 2026_06_16_194921, 41 files, 0 failures, 465 s); each band
is a 2-universe band → `pair_covariance`, so
`cov_RPA = pair(HighQ2⁺,HighQ2⁻) + pair(LowQ2⁺,LowQ2⁻)`.

The error-band construction (`RPAWeight._weights_one_z`) is a port of
weightRPA.cxx:55-235, **parity-tested** by hand-replicating the band math from
the `hrelratio`/`hnonrelratio`/`hQ2rel`/`hQ2nonrel` histograms; the CV path is
verified unchanged to 1e-12.

## Result (playlist 1A)

| metric | value |
|---|---|
| RPA fractional uncertainty / cell (median) | **0.097 %** |
| HighQ2 (enhancement) / LowQ2 (suppression) | 0.043 % / 0.037 % |
| p84 / max | 0.22 % / 1.03 % |
| **integrated σ across the 4 universes** | 1.000 / 1.000 / 0.9999 / 1.0001 |
| cov_RPA symmetric / PSD | yes / yes |

**RPA is negligible for this inclusive measurement** (~0.1 %/cell, integrated σ
unchanged to ≤0.01 %). The reason: it touches only true-QE on Z≥6 — a subset of
the inclusive sample — and enters solely through the MC efficiency/migration
ratios, which are data-normalized and largely cancel. (By contrast it is a
*leading* systematic for the QE-exclusive companion analysis, where the QE
content and its Q² shape are the measurement.)

## Effect on the interaction-model category

| interaction-model band | frac/cell (1A) |
|---|---|
| GENIE (56 universes) | 1.19 % |
| 2p2h (3 universes) | 0.31 % |
| **RPA (4 universes)** | **0.10 %** |

The systematic-only total is unchanged at 5.30 % (RPA adds in quadrature
negligibly). The remaining interaction-model piece is **GenieRvx1pi** (the
non-resonant single-π normalization, confirmed a separate band in the
`GetStandardSystematics` dispatch) — the likely source of the residual gap to
the Fig 8 "Models" curve (~2–3 %).

## Reproduce

```bash
pixi run python make_rpa_universes.py --workers 8 --playlist minervame1A
pixi run python assemble_total.py --ingredients <ing> --xsec <xsec> --muon <cov_muon> \
    --stat-cov <stat> --genie <genie> --twop2h <twop2h> --rpa <rpa>
```
