# POT ledger — all 12 ME FHC playlists

`pot_ledger.py --workers 12 --label all12_playlists` (2026-06-12): streamed the
Meta tree of **all 2,307 files** (1,818 data + 489 StandardMC), **0 failures,
43.6 s**. RunLog `~/log/ndp-minerva-xsec/2026_06_12_211937.log`; full per-file
ledger in `results/<ts>__pot_ledger/ledger.json` (untracked).

## Headline

**Σ POT_Used (data, 1A–1P) = 1.0574×10²¹ = 99.66 % of the paper's 10.61×10²⁰.**
The open-data release carries essentially exactly the paper's exposure — use
all 12 playlists; no run-range subset hunting is needed. The long-standing
"10.61 vs 12.13/11.12" accounting puzzle is closed:

- the 12.13×10²⁰ figure (early planning doc) was simply wrong;
- the getdata-page per-playlist values are coarse **delivered-POT-style**
  numbers: their sum (11.12×10²⁰) matches our Σ POT_Total (11.079×10²⁰) to
  0.4 %;
- **minervame1M is the single genuine anomaly**: preserved files hold
  POT_Used 1.581 / POT_Total 1.714 vs page 2.1 (×10²⁰) — ~25 % of the period
  is not in the preserved set (consistent with NOT-PRESERVED runs and 1M's
  missing edge runs 19185/19484). All other playlists agree with their page
  value within ±3 % (POT_Used, allowing the page's 2-s.f. rounding).

## Per playlist

| playlist | data files | data POT_Used | data POT_Total | page (2 s.f.) | used/page | mc files | mc POT_Used | mc/data |
|---|---|---|---|---|---|---|---|---|
| 1A | 253 | 8.9689e19 | 9.6028e19 | 0.90e20 | 0.997 | 41 | 4.0717e20 | 4.54 |
| 1B | 47 | 1.8674e19 | 1.8792e19 | 0.19e20 | 0.983 | 11 | 1.0937e20 | 5.86 |
| 1C | 102 | 4.2944e19 | 4.9566e19 | 0.43e20 | 0.999 | 21 | 2.0874e20 | 4.86 |
| 1D | 283 | 1.4404e20 | 1.4506e20 | 1.4e20 | 1.029 | 61 | 6.0740e20 | 4.22 |
| 1E | 219 | 1.0296e20 | 1.0507e20 | 1.0e20 | 1.030 | 51 | 5.0861e20 | 4.94 |
| 1F | 260 | 1.6686e20 | 1.7225e20 | 1.7e20 | 0.982 | 71 | 7.0712e20 | 4.24 |
| 1G | 233 | 1.3757e20 | 1.4621e20 | 1.4e20 | 0.983 | 61 | 5.9543e20 | 4.33 |
| 1L | 15 | 1.3371e19 | 1.3652e19 | 0.13e20 | 1.029 | 6 | 5.8228e19 | 4.35 |
| **1M** | 202 | **1.5813e20** | 1.7143e20 | **2.1e20** | **0.753** | 98 | 8.9789e20 | 5.68 |
| 1N | 138 | 1.0657e20 | 1.1251e20 | 1.1e20 | 0.969 | 31 | 5.1211e20 | 4.81 |
| 1O | 26 | 2.9797e19 | 2.9913e19 | 0.30e20 | 0.993 | 16 | 1.5824e20 | 5.31 |
| 1P | 40 | 4.6758e19 | 4.7414e19 | 0.47e20 | 0.995 | 21 | 2.0814e20 | 4.45 |
| **total** | **1818** | **1.0574e21** | 1.1079e21 | (11.12e20) | — | **489** | **4.9784e21** | **4.71** |

## Consequences for the analysis

1. **Normalization input settled** (pending user confirmation): data POT for
   the full-dataset extraction = 1.0574×10²¹ (POT_Used over all 12 playlists);
   matches the paper's normalization to 0.34 %, far below the 3.9 % flux
   normalization uncertainty.
2. MC exposure totals 4.978×10²¹ (4.71× data; per-playlist 4.2–5.9) — fixes
   the per-playlist background/efficiency scaling for the multi-playlist pass.
3. Per-playlist POT pairs are now recorded per file in `ledger.json` — the
   future skim/ingredients passes can carry exact POT without re-reading Meta.

## Reproduce

```bash
pixi run python pot_ledger.py --workers 12 --label all12_playlists
```
