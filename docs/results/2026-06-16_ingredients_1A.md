# E1 — six cross-section ingredients, playlist 1A (weighted CV)

`make_ingredients.py --weights cv --playlist minervame1A` (RunLog
2026_06_16_142140): streams all 253 data + 41 MC files, both the reco tree
(data hist, background, signal migration) and the **MC Truth tree** (efficiency
denominator — the new loop), into count-conserving slot space
(`xsec/binning.py`, N_SLOTS=288). Output: `results/<ts>__make_ingredients/
ingredients.npz` (data_reco, bkg, migration[reco,true], eff_num, eff_denom,
POT, slot maps).

## Gates — all pass

| gate | value | reference |
|---|---|---|
| data selected, in grid | 357,547 | = weighted-1A run ✓ |
| MC signal selected (total) | 1,786,201 | = weighted-1A run ✓ |
| migration total (count-conserving) | 1,453,610 | = weighted sumw_signal; 0 dropped ✓ |
| eff_denom (unweighted total) | 2,684,465 | = 41 × 65,712 ✓ |
| **average efficiency (weighted)** | **0.657** | paper ≈ 0.64 ✓ |
| single-file anchor (unweighted) | eff_denom **65,712** | frozen exploration anchor, exact ✓ |

The efficiency asymmetry is the tutorial Cutter convention (verified):
numerator = reco-selected signal (no truth phase-space cut), free as the
migration's reco-axis column-sum; denominator = signal ∧ phase-space from the
Truth tree. Reco fills carry the full reco CV weight (incl. MINOS-eff); the
Truth denominator carries the truth-only weight — matching Model semantics.

## Notes

- Streaming-only; nothing saved beyond the histograms. The `ingredients.npz`
  is the sole input to E2 (unfolding) — no re-streaming downstream.
- Runtime caveat: the full run wall time was dominated by a single ~14.7 h
  xrootd door stall between files 150–175 (cadence was ~0.6 s/file before and
  after; profiled per-file weighted-truth work is ~7 s). A per-read xrootd
  timeout + the existing retry would re-dispatch a hung door instead of
  blocking the pool — a robustness follow-up, not a correctness issue.

## Reproduce

```bash
pixi run python make_ingredients.py --weights cv --playlist minervame1A
```
