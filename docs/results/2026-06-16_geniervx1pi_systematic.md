# GenieRvx1pi — non-resonant single-π normalization (interaction model)

The MnvTune non-res-π reduction (CV weight **0.43** on tagged events) carries an
uncertainty, `GenieRvx1pi`. The dispatch `GetGenieRvx1piSystematicsMap`
registers **two bands** — "Rvn1pi" and "Rvp1pi", each a ±1σ pair. With the
non-res-π reweight applied (our MnvTune v1 case), both bands **collapse to the
same flat shift** (`GenieRvx1piUniverse::GetWeightRatioToCV`, reweight-on
branch, GenieSystematics.cxx:674-687):

  ratio to CV = 1 ± kNonResPiWeightShift / kNonResPiWeight = 1 ± 0.04/0.43 ≈ **±9.3 %**

on the `IsNonResPi`-tagged events (`Rv{n,p}1pi[2] < 1`, identical to the CV tag)
and 1 elsewhere. So the universe weight on a tagged event is 0.43 ± 0.04.

`make_geniervx1pi_universes.py` runs the ±1σ vertical pair (RunLog
2026_06_16_214932, 41 files, 0 failures, 432 s). Because the two MAT bands are
identical under the reweight, the group covariance is their **sum = 2 ×
pair_covariance(σ⁺, σ⁻)** — `assemble_total.py --geniervx1pi` applies the 2×
and reports both the 1- and 2-band magnitudes.

## Result (playlist 1A)

| metric | value |
|---|---|
| GenieRvx1pi **2-band** (folded into total) | **0.026 %/cell** (median) |
| GenieRvx1pi **1-band** (single pair) | 0.018 %/cell |
| 2-band / 1-band | √2 (= the two identical bands) |
| p84 / p98 / max | 0.079 % / 0.38 % / 0.63 % |
| cells > 0.5 % / > 1 % | 2 / 0 (of 205) |
| symmetric / PSD | yes / yes |

**Negligible for the inclusive measurement.** It is ~0 in most cells (median
0.026 %) and reaches only ~0.6 % in the high-p_T (DIS) and high-p_∥ corners
where non-resonant pion production lives. The flat ±9.3 % on the tagged subset
dilutes and largely cancels in the data-driven extraction (conditional
migration + eff num/denom). The **2× double-counting subtlety is moot here** —
the band is negligible whether counted as 1× or 2×.

## Interaction-model category — now complete

| band | frac/cell (1A) |
|---|---|
| GENIE (56 knob universes) | 1.19 % |
| 2p2h (3 universes) | 0.31 % |
| RPA (4 universes) | 0.10 % |
| **GenieRvx1pi (non-res π)** | **0.03 %** |

All four interaction-model bands of `GetStandardSystematics` are now built, and
the category is **GENIE-dominated at ~1.2 %/cell** — consistent with the paper's
text ("the interaction model is not dominant in any (p_T, p_∥) bin"). This sits
below the Fig 8 "Models" curve's apparent ~2–3 %; the residual is most likely
the **"new" GENIE systematics** (the deuterium-fit MaRES⊗NormCCRES covariance
band and the FaCCQE z-expansion, `GenieMaNormResCovUniverse`/`GenieFaCCQEUniverse`)
which replace the simple knobs in newer MAT — a candidate follow-up, distinct
from the `GenieRvx1pi` band done here.

## Effect on Cov_total (1A)
Systematic-only unchanged at **5.30 %** (GenieRvx1pi adds in quadrature
negligibly). The remaining groups are geant4/response (Hadronic Response),
flux shape, and a normalization band.

## Reproduce
```bash
pixi run python make_geniervx1pi_universes.py --workers 8 --playlist minervame1A
pixi run python assemble_total.py --ingredients <ing> --xsec <xsec> --muon <cov_muon> \
    --stat-cov <stat> --genie <genie> --twop2h <twop2h> --rpa <rpa> --geniervx1pi <grvx>
```
