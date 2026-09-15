# Selection efficiency maps and background of the MINERvA 1μ1p selection (full ME FHC official MC)

**Status:** rendered 2026-09-15 from `runs/2026-09-15_efficiency_minerva_me_ccqelike_1mu1p` (`ndp efficiency run`), which reads the binned responses built by
`ndp surrogate build --channel minerva_me_ccqelike_1mu1p --measurement all` over the 12 FHC playlists
(`surrogates/minerva_me_ccqelike_1mu1p/<grid>/binned_FHC_1A-1P/`, every closure exact). Per-grid tables live in that run directory:
`eff_<grid>.json` / `.npz` (edges, efficiency, binomial error, denominator, numerator), `background_<grid>.json`
(total, by category, feed-in, at the data POT) and `ansatz_closure.json` (the per-bin closure of the factorised maps).

**How to use this on a truth sample.** For a weight in the bins of one released variable, take that grid's per-bin ε
below (exact for its own binning). For a per-event weight from the two maps,

```
python -m ndp efficiency apply --channel minerva_me_ccqelike_1mu1p --run runs/2026-09-15_efficiency_minerva_me_ccqelike_1mu1p --sample <truth.npz> --out weights.npz
```

writes w = ε_μ(p_μ, cos θ_μ) · ε_p(p_p, cos θ_p) / ⟨ε⟩ per event (0 outside the maps or without a leading proton in
the window). The closure table says how well that factorised weight reproduces the actually selected signal on each
grid: within a few per cent on the muon and proton kinematics, but far off on the transverse-imbalance variables,
where the selection efficiency depends on the muon–proton correlation. For those, fold through the full response
(`python -m ndp run <model> --channel minerva_me_ccqelike_1mu1p --measurement all --modes folded`).

---


Official MC 4.978e+21 POT: 2734225 truth signal events in the fiducial volume, 891177 of them selected: mean efficiency 0.3259. Backgrounds are scaled to the data POT 1.057e+21. Efficiency = selected truth signal / truth signal in the fiducial volume per true cell (the paper's definition); error = binomial.

## Efficiency per grid

| measurement | grid | <eff> | den | background at data POT |
|---|---|---|---|---|
| muon_p_costheta | 8 × 8 | 0.3259 | 2734225 | 202110 |
| proton_p_costheta | 5 × 8 | 0.3259 | 2734227 | 166633 |
| muon_costheta | 8 | 0.3259 | 2734225 | 202110 |
| proton_costheta | 8 | 0.3259 | 2734227 | 198410 |
| alpha | 5 | 0.3259 | 2734227 | 201982 |
| dpt | 6 | 0.3259 | 2734207 | 201827 |
| dpt_fine | 9 | 0.3259 | 2734207 | 201827 |
| dptx | 8 | 0.3259 | 2734227 | 201982 |
| dpty | 12 | 0.3259 | 2734227 | 201982 |
| muon_p | 8 | 0.3259 | 2734227 | 202110 |
| muon_pt | 10 | 0.3259 | 2734227 | 202105 |
| muon_theta | 14 | 0.3259 | 2734227 | 202110 |
| phi | 6 | 0.3259 | 2734227 | 201982 |
| pl | 8 | 0.3259 | 2734227 | 201947 |
| pn | 10 | 0.3259 | 2734227 | 201982 |
| proton_p | 5 | 0.3259 | 2734227 | 170230 |
| proton_pt | 7 | 0.3259 | 2734227 | 201926 |
| proton_theta | 10 | 0.3259 | 2734227 | 198410 |

### muon_p_costheta: ε(muon p, cos theta_mu)

| muon p \ cos theta_mu | [0.9563, 0.9618) | [0.9618, 0.9672) | [0.9672, 0.9727) | [0.9727, 0.9782) | [0.9782, 0.9836) | [0.9836, 0.9891) | [0.9891, 0.9945) | [0.9945, 1) |
|---|---|---|---|---|---|---|---|---|
| [2, 3) | 0.101 ± 0.002 | 0.147 ± 0.002 | 0.183 ± 0.002 | 0.222 ± 0.002 | 0.259 ± 0.002 | 0.292 ± 0.002 | 0.342 ± 0.003 | 0.408 ± 0.003 |
| [3, 4) | 0.096 ± 0.002 | 0.135 ± 0.002 | 0.172 ± 0.002 | 0.221 ± 0.002 | 0.273 ± 0.002 | 0.318 ± 0.002 | 0.354 ± 0.002 | 0.416 ± 0.002 |
| [4, 5) | 0.073 ± 0.003 | 0.108 ± 0.003 | 0.149 ± 0.003 | 0.183 ± 0.002 | 0.242 ± 0.002 | 0.305 ± 0.001 | 0.361 ± 0.001 | 0.402 ± 0.001 |
| [5, 6) | 0.050 ± 0.005 | 0.087 ± 0.005 | 0.119 ± 0.004 | 0.163 ± 0.003 | 0.207 ± 0.002 | 0.276 ± 0.001 | 0.356 ± 0.001 | 0.403 ± 0.001 |
| [6, 7.5) | 0.037 ± 0.012 | 0.056 ± 0.009 | 0.067 ± 0.007 | 0.123 ± 0.006 | 0.176 ± 0.004 | 0.246 ± 0.002 | 0.338 ± 0.001 | 0.407 ± 0.001 |
| [7.5, 10) | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.041 ± 0.018 | 0.061 ± 0.014 | 0.117 ± 0.010 | 0.189 ± 0.006 | 0.301 ± 0.002 | 0.403 ± 0.001 |
| [10, 14) | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.069 ± 0.047 | 0.034 ± 0.019 | 0.083 ± 0.014 | 0.186 ± 0.008 | 0.389 ± 0.002 |
| [14, 20) | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.061 ± 0.034 | 0.064 ± 0.013 | 0.351 ± 0.003 |

### proton_p_costheta: ε(leading proton p, cos theta_p)

| leading proton p \ cos theta_p | [0.342, 0.4243) | [0.4243, 0.5065) | [0.5065, 0.5888) | [0.5888, 0.671) | [0.671, 0.7533) | [0.7533, 0.8355) | [0.8355, 0.9178) | [0.9178, 1) |
|---|---|---|---|---|---|---|---|---|
| [0.5, 0.625) | 0.114 ± 0.001 | 0.217 ± 0.001 | 0.296 ± 0.001 | 0.341 ± 0.002 | 0.350 ± 0.002 | 0.370 ± 0.002 | 0.377 ± 0.002 | 0.382 ± 0.002 |
| [0.625, 0.75) | 0.307 ± 0.001 | 0.377 ± 0.001 | 0.413 ± 0.002 | 0.432 ± 0.002 | 0.441 ± 0.002 | 0.452 ± 0.002 | 0.458 ± 0.002 | 0.428 ± 0.002 |
| [0.75, 0.875) | 0.283 ± 0.002 | 0.352 ± 0.001 | 0.386 ± 0.001 | 0.404 ± 0.002 | 0.408 ± 0.002 | 0.425 ± 0.002 | 0.427 ± 0.002 | 0.405 ± 0.002 |
| [0.875, 1) | 0.219 ± 0.002 | 0.255 ± 0.002 | 0.278 ± 0.001 | 0.290 ± 0.002 | 0.309 ± 0.002 | 0.332 ± 0.002 | 0.351 ± 0.002 | 0.351 ± 0.002 |
| [1, 1.1) | 0.153 ± 0.003 | 0.165 ± 0.002 | 0.178 ± 0.002 | 0.195 ± 0.002 | 0.215 ± 0.002 | 0.246 ± 0.003 | 0.274 ± 0.003 | 0.282 ± 0.003 |

## Efficiency per bin of every one-dimensional grid

These are the numbers to weight a truth sample with when the weight is wanted in the bins of one released variable: ε of a true bin is the fraction of that bin's fiducial signal events the selection keeps, so a prediction in that bin is (truth events in the bin) × ε. They are exact for their own binning, unlike the factorised map product below.

**muon_costheta** — ε vs true cos theta_mu

| bin | [0.9563, 0.9618) | [0.9618, 0.9672) | [0.9672, 0.9727) | [0.9727, 0.9782) | [0.9782, 0.9836) | [0.9836, 0.9891) | [0.9891, 0.9945) | [0.9945, 1) |
|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.094 ± 0.001 | 0.134 ± 0.001 | 0.167 ± 0.001 | 0.205 ± 0.001 | 0.248 ± 0.001 | 0.293 ± 0.001 | 0.349 ± 0.001 | 0.403 ± 0.000 |
| signal in bin | 56368 | 76434 | 105489 | 152103 | 232505 | 391528 | 714480 | 1005318 |

**proton_costheta** — ε vs true cos theta_p

| bin | [0.342, 0.4243) | [0.4243, 0.5065) | [0.5065, 0.5888) | [0.5888, 0.671) | [0.671, 0.7533) | [0.7533, 0.8355) | [0.8355, 0.9178) | [0.9178, 1) |
|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.214 ± 0.001 | 0.294 ± 0.001 | 0.325 ± 0.001 | 0.342 ± 0.001 | 0.359 ± 0.001 | 0.382 ± 0.001 | 0.392 ± 0.001 | 0.383 ± 0.001 |
| signal in bin | 392978 | 484398 | 460094 | 370549 | 279384 | 254321 | 250087 | 242416 |

**alpha** — ε vs true δα_T [deg]

| bin | [0, 40) | [40, 80) | [80, 120) | [120, 150) | [150, 180) |
|---|---|---|---|---|---|
| ε ± δ | 0.345 ± 0.001 | 0.339 ± 0.001 | 0.339 ± 0.001 | 0.328 ± 0.001 | 0.306 ± 0.000 |
| signal in bin | 382941 | 313246 | 463057 | 620719 | 954264 |

**dpt** — ε vs true δp_T (coarse) [GeV/c]

| bin | [0, 0.2) | [0.2, 0.4) | [0.4, 0.6) | [0.6, 1) | [1, 1.75) | [1.75, 3) |
|---|---|---|---|---|---|---|
| ε ± δ | 0.327 ± 0.000 | 0.364 ± 0.001 | 0.350 ± 0.001 | 0.307 ± 0.001 | 0.199 ± 0.001 | 0.083 ± 0.003 |
| signal in bin | 1148181 | 560067 | 391016 | 460173 | 168139 | 6631 |

**dpt_fine** — ε vs true δp_T [GeV/c]

| bin | [0, 0.1) | [0.1, 0.2) | [0.2, 0.3) | [0.3, 0.4) | [0.4, 0.5) | [0.5, 0.6) | [0.6, 0.8) | [0.8, 1.5) | [1.5, 3) |
|---|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.319 ± 0.001 | 0.332 ± 0.001 | 0.361 ± 0.001 | 0.369 ± 0.001 | 0.357 ± 0.001 | 0.342 ± 0.001 | 0.323 ± 0.001 | 0.247 ± 0.001 | 0.112 ± 0.002 |
| signal in bin | 501555 | 646626 | 315400 | 244667 | 210522 | 180494 | 282051 | 330056 | 22836 |

**dptx** — ε vs true δp_Tx [GeV/c]

| bin | [-2, -0.7) | [-0.7, -0.4) | [-0.4, -0.2) | [-0.2, 0) | [0, 0.2) | [0.2, 0.4) | [0.4, 0.7) | [0.7, 2) |
|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.217 ± 0.003 | 0.291 ± 0.001 | 0.337 ± 0.001 | 0.331 ± 0.001 | 0.330 ± 0.001 | 0.339 ± 0.001 | 0.293 ± 0.001 | 0.220 ± 0.003 |
| signal in bin | 16835 | 174748 | 295119 | 880215 | 879423 | 296083 | 174604 | 17200 |

**dpty** — ε vs true δp_Ty [GeV/c]

| bin | [-6.5, -3.5) | [-3.5, -2) | [-2, -1.1) | [-1.1, -0.8) | [-0.8, -0.6) | [-0.6, -0.4) | [-0.4, -0.2) | [-0.2, 0) | [0, 0.2) | [0.2, 0.4) | [0.4, 1.3) | [1.3, 3) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.000 ± 0.000 | 0.061 ± 0.006 | 0.167 ± 0.001 | 0.261 ± 0.001 | 0.312 ± 0.001 | 0.341 ± 0.001 | 0.354 ± 0.001 | 0.327 ± 0.001 | 0.341 ± 0.001 | 0.357 ± 0.002 | 0.313 ± 0.004 | 0.000 ± 0.000 |
| signal in bin | 2 | 1500 | 91541 | 168777 | 195354 | 267662 | 380902 | 838840 | 690570 | 83098 | 15981 | 0 |

**muon_p** — ε vs true muon p [GeV/c]

| bin | [2, 3) | [3, 4) | [4, 5) | [5, 6) | [6, 7.5) | [7.5, 10) | [10, 14) | [14, 20) |
|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.236 ± 0.001 | 0.284 ± 0.001 | 0.317 ± 0.001 | 0.344 ± 0.001 | 0.368 ± 0.001 | 0.382 ± 0.001 | 0.377 ± 0.002 | 0.346 ± 0.003 |
| signal in bin | 264185 | 470127 | 582864 | 578697 | 538938 | 223608 | 51491 | 24317 |

**muon_pt** — ε vs true muon p_T [GeV/c]

| bin | [0, 0.1) | [0.1, 0.2) | [0.2, 0.3) | [0.3, 0.45) | [0.45, 0.6) | [0.6, 0.75) | [0.75, 0.9) | [0.9, 1.25) | [1.25, 2.5) | [2.5, 5) |
|---|---|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.456 ± 0.004 | 0.448 ± 0.002 | 0.430 ± 0.001 | 0.379 ± 0.001 | 0.357 ± 0.001 | 0.326 ± 0.001 | 0.279 ± 0.001 | 0.210 ± 0.001 | 0.113 ± 0.002 | 0.000 ± 0.000 |
| signal in bin | 19758 | 66858 | 126032 | 374276 | 641671 | 699958 | 476614 | 294660 | 34334 | 66 |

**muon_theta** — ε vs true muon θ [deg]

| bin | [0, 1) | [1, 2) | [2, 3) | [3, 4) | [4, 5) | [5, 6) | [6, 7) | [7, 8) | [8, 9) | [9, 10) | [10, 11) | [11, 12) | [12, 14) | [14, 17) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.459 ± 0.004 | 0.445 ± 0.002 | 0.427 ± 0.001 | 0.407 ± 0.001 | 0.395 ± 0.001 | 0.384 ± 0.001 | 0.364 ± 0.001 | 0.342 ± 0.001 | 0.317 ± 0.001 | 0.290 ± 0.001 | 0.265 ± 0.001 | 0.239 ± 0.001 | 0.198 ± 0.001 | 0.129 ± 0.001 |
| signal in bin | 20058 | 69259 | 132740 | 205331 | 272027 | 308578 | 308838 | 281880 | 242892 | 201530 | 165615 | 135230 | 202842 | 187407 |

**phi** — ε vs true φ_T [deg]

| bin | [0, 10) | [10, 30) | [30, 50) | [50, 80) | [80, 120) | [120, 180) |
|---|---|---|---|---|---|---|
| ε ± δ | 0.315 ± 0.000 | 0.338 ± 0.001 | 0.341 ± 0.001 | 0.333 ± 0.001 | 0.324 ± 0.001 | 0.316 ± 0.001 |
| signal in bin | 1056984 | 624950 | 291184 | 276419 | 235651 | 249039 |

**pl** — ε vs true δp_L [GeV/c]

| bin | [-1, -0.25) | [-0.25, 0) | [0, 0.1) | [0.1, 0.2) | [0.2, 0.3) | [0.3, 0.4) | [0.4, 0.6) | [0.6, 1) |
|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.060 ± 0.003 | 0.220 ± 0.001 | 0.299 ± 0.001 | 0.366 ± 0.001 | 0.399 ± 0.001 | 0.436 ± 0.001 | 0.407 ± 0.002 | 0.000 ± 0.000 |
| signal in bin | 5562 | 603097 | 818442 | 565050 | 404135 | 260848 | 77093 | 0 |

**pn** — ε vs true p_n [GeV/c]

| bin | [0, 0.1) | [0.1, 0.2) | [0.2, 0.3) | [0.3, 0.4) | [0.4, 0.5) | [0.5, 0.6) | [0.6, 0.7) | [0.7, 1) | [1, 2) | [2, 6) |
|---|---|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.290 ± 0.001 | 0.318 ± 0.001 | 0.332 ± 0.001 | 0.378 ± 0.001 | 0.386 ± 0.001 | 0.367 ± 0.001 | 0.350 ± 0.001 | 0.309 ± 0.001 | 0.202 ± 0.001 | 0.059 ± 0.005 |
| signal in bin | 307164 | 635919 | 308901 | 257292 | 271133 | 227460 | 185280 | 354562 | 184595 | 1921 |

**proton_p** — ε vs true leading proton p [GeV/c]

| bin | [0.5, 0.625) | [0.625, 0.75) | [0.75, 0.875) | [0.875, 1) | [1, 1.1) |
|---|---|---|---|---|---|
| ε ± δ | 0.278 ± 0.001 | 0.402 ± 0.001 | 0.380 ± 0.001 | 0.292 ± 0.001 | 0.207 ± 0.001 |
| signal in bin | 725712 | 687646 | 576240 | 461475 | 283154 |

**proton_pt** — ε vs true leading proton p_T [GeV/c]

| bin | [0, 0.15) | [0.15, 0.3) | [0.3, 0.45) | [0.45, 0.6) | [0.6, 0.75) | [0.75, 0.9) | [0.9, 1.1) |
|---|---|---|---|---|---|---|---|
| ε ± δ | 0.374 ± 0.002 | 0.392 ± 0.001 | 0.371 ± 0.001 | 0.319 ± 0.001 | 0.345 ± 0.001 | 0.238 ± 0.001 | 0.156 ± 0.001 |
| signal in bin | 67174 | 218611 | 455727 | 867489 | 675833 | 385243 | 64150 |

**proton_theta** — ε vs true leading proton θ [deg]

| bin | [0, 10) | [10, 16.25) | [16.25, 22.5) | [22.5, 28.75) | [28.75, 35) | [35, 41.25) | [41.25, 47.5) | [47.5, 53.75) | [53.75, 60) | [60, 70) |
|---|---|---|---|---|---|---|---|---|---|---|
| ε ± δ | 0.363 ± 0.002 | 0.380 ± 0.002 | 0.391 ± 0.001 | 0.395 ± 0.001 | 0.389 ± 0.001 | 0.380 ± 0.001 | 0.359 ± 0.001 | 0.343 ± 0.001 | 0.324 ± 0.001 | 0.256 ± 0.000 |
| signal in bin | 44158 | 72598 | 107537 | 142213 | 176548 | 208123 | 257000 | 375630 | 512120 | 838300 |


## Factorised ansatz w = ε_μ(p_μ, cos θ_μ) · ε_p(p_p, cos θ_p) / ⟨ε⟩ — closure on the MC

⟨ε⟩ = 0.3259; total ansatz-selected / actually selected = 1.0032.

| grid | ansatz / actual (total) | max |rel. dev.| over bins with > 50 events |
|---|---|---|
| muon_p_costheta | 1.0032 | 0.098 |
| proton_p_costheta | 1.0032 | 0.133 |
| muon_costheta | 1.0032 | 0.022 |
| proton_costheta | 1.0032 | 0.012 |
| alpha | 1.0032 | 0.039 |
| dpt | 1.0032 | 0.458 |
| dpt_fine | 1.0032 | 0.398 |
| dptx | 1.0032 | 0.124 |
| dpty | 1.0032 | 0.507 |
| muon_p | 1.0032 | 0.043 |
| muon_pt | 1.0032 | 0.237 |
| muon_theta | 1.0032 | 0.054 |
| phi | 1.0032 | 0.164 |
| pl | 1.0032 | 0.458 |
| pn | 1.0032 | 0.582 |
| proton_p | 1.0032 | 0.088 |
| proton_pt | 1.0032 | 0.061 |
| proton_theta | 1.0032 | 0.070 |

The ansatz is exact by construction for the two map grids' totals; the per-grid deviations are the error made by weighting a truth sample with the maps instead of folding it through the full response.

## Background composition (selected non-signal MC at the data POT)

| grid | total | bkg 1 pi+- | bkg 1 pi0 | bkg multi-pi | bkg other (no pion) | feed-in |
|---|---|---|---|---|---|---|
| muon_p_costheta | 202110 | 94927 | 45622 | 25375 | 36185 | 0 |
| proton_p_costheta | 166633 | 83791 | 40141 | 22151 | 20551 | 0 |
| muon_costheta | 202110 | 94927 | 45622 | 25375 | 36185 | 0 |
| proton_costheta | 198410 | 94477 | 45456 | 25320 | 33157 | 0 |
| alpha | 201982 | 94892 | 45609 | 25371 | 36111 | 0 |
| dpt | 201827 | 94840 | 45570 | 25334 | 36082 | 0 |
| dpt_fine | 201827 | 94840 | 45570 | 25334 | 36082 | 0 |
| dptx | 201982 | 94892 | 45609 | 25371 | 36111 | 0 |
| dpty | 201982 | 94892 | 45609 | 25371 | 36111 | 0 |
| muon_p | 202110 | 94927 | 45622 | 25375 | 36185 | 0 |
| muon_pt | 202105 | 94925 | 45622 | 25374 | 36184 | 0 |
| muon_theta | 202110 | 94927 | 45622 | 25375 | 36185 | 0 |
| phi | 201982 | 94892 | 45609 | 25371 | 36111 | 0 |
| pl | 201947 | 94879 | 45600 | 25366 | 36103 | 0 |
| pn | 201982 | 94891 | 45609 | 25371 | 36110 | 0 |
| proton_p | 170230 | 84202 | 40297 | 22200 | 23531 | 0 |
| proton_pt | 201926 | 94870 | 45603 | 25367 | 36085 | 0 |
| proton_theta | 198410 | 94477 | 45456 | 25320 | 33157 | 0 |

## Figures

![eff_ansatz_alpha](figs/eff_ansatz_alpha.png)
![eff_ansatz_dpt](figs/eff_ansatz_dpt.png)
![eff_ansatz_dpt_fine](figs/eff_ansatz_dpt_fine.png)
![eff_ansatz_dptx](figs/eff_ansatz_dptx.png)
![eff_ansatz_dpty](figs/eff_ansatz_dpty.png)
![eff_ansatz_muon_costheta](figs/eff_ansatz_muon_costheta.png)
![eff_ansatz_muon_p](figs/eff_ansatz_muon_p.png)
![eff_ansatz_muon_pt](figs/eff_ansatz_muon_pt.png)
![eff_ansatz_muon_theta](figs/eff_ansatz_muon_theta.png)
![eff_ansatz_phi](figs/eff_ansatz_phi.png)
![eff_ansatz_pl](figs/eff_ansatz_pl.png)
![eff_ansatz_pn](figs/eff_ansatz_pn.png)
![eff_ansatz_proton_costheta](figs/eff_ansatz_proton_costheta.png)
![eff_ansatz_proton_p](figs/eff_ansatz_proton_p.png)
![eff_ansatz_proton_pt](figs/eff_ansatz_proton_pt.png)
![eff_ansatz_proton_theta](figs/eff_ansatz_proton_theta.png)
![eff_eff_map_muon_p_costheta](figs/eff_eff_map_muon_p_costheta.png)
![eff_eff_map_proton_p_costheta](figs/eff_eff_map_proton_p_costheta.png)
![eff_eff_proj_alpha](figs/eff_eff_proj_alpha.png)
![eff_eff_proj_dpt](figs/eff_eff_proj_dpt.png)
![eff_eff_proj_dpt_fine](figs/eff_eff_proj_dpt_fine.png)
![eff_eff_proj_dptx](figs/eff_eff_proj_dptx.png)
![eff_eff_proj_dpty](figs/eff_eff_proj_dpty.png)
![eff_eff_proj_muon_costheta](figs/eff_eff_proj_muon_costheta.png)
![eff_eff_proj_muon_p](figs/eff_eff_proj_muon_p.png)
![eff_eff_proj_muon_p_costheta](figs/eff_eff_proj_muon_p_costheta.png)
![eff_eff_proj_muon_pt](figs/eff_eff_proj_muon_pt.png)
![eff_eff_proj_muon_theta](figs/eff_eff_proj_muon_theta.png)
![eff_eff_proj_phi](figs/eff_eff_proj_phi.png)
![eff_eff_proj_pl](figs/eff_eff_proj_pl.png)
![eff_eff_proj_pn](figs/eff_eff_proj_pn.png)
![eff_eff_proj_proton_costheta](figs/eff_eff_proj_proton_costheta.png)
![eff_eff_proj_proton_p](figs/eff_eff_proj_proton_p.png)
![eff_eff_proj_proton_p_costheta](figs/eff_eff_proj_proton_p_costheta.png)
![eff_eff_proj_proton_pt](figs/eff_eff_proj_proton_pt.png)
![eff_eff_proj_proton_theta](figs/eff_eff_proj_proton_theta.png)
