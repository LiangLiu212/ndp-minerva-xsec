# GiBUU 1μ1p signal on the MINERvA ME FHC selection: efficiency, background and data overlay

**Status:** rendered from `runs/2026-09-15_efficiency_minerva_me_ccqelike_1mu1p` (efficiency maps, background, ansatz closure) and `runs/2026-09-16_gibuu_2025_me_fhc_c12__minerva_me_ccqelike_1mu1p__all_3` (data vs GiBUU signal + MC background on 18 grids), manifest git `2b81262`. Numbers are quoted from those directories; this file was rendered by `report/make_overlay_report.py`.

**Model:** GiBUU 2025 numu CC on 12C under the NuMI ME FHC flux, generated as fixed-energy integratedSigma runs on a flux-weighted 0.5 GeV energy grid (2-60 GeV, 116 points, 156 grid jobs, ~1.07M events with near-uniform weights); in-medium resonance widths on, 2p2h and two-pion background on.

---

## 1. Inputs

| item | value |
|---|---|
| data | ME FHC playlists 1A–1P, 1.057 × 10²¹ POT_Used, 329 653 selected events |
| official MC (efficiency and background) | 4.978 × 10²¹ POT, 2 734 225 truth signal events in the fiducial volume, 891 177 selected (⟨ε⟩ = 0.3259) |
| GiBUU sample | GiBUU Release 2025 patch 5, 979 395 events from 153 job(s), fingerprint `1f8bcb950be1dcd0` |
| GiBUU σ_CC (flux-averaged, per nucleon) | 4.1245e-38 cm² |
| normalisation | N_true = sigma_cell * N_nuc(3.23e+30) * Phi(6.32e-08) * POT(1.057e+21) |

## 1b. The GiBUU sample

115 energy points from 2.75 to 59.75 GeV covering 0.9460 of the 0 to 100 GeV flux, 153 grid jobs, 979 395 events. Flux-averaged σ$_{CC}$ over the covered flux: 4.1245e-38 cm²/nucleon. **1 planned points have no job**, so the predicted rate is biased low by their share of the flux.

![sigma of energy](figs/gibuu_sigma_of_energy.png)
*σ$_{CC}$(E) of every point against the GENIE spline, with each point's flux weight and event count.*

## 2. Selection efficiency (official MC)

ε = selected truth signal / truth signal in the fiducial volume, per true cell (the paper's definition), binomial errors. The two maps are what a truth sample can be weighted with; the 1D grids are their projections and the released grids' efficiencies.

### ε(muon p, cos theta_mu)

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

![muon_p_costheta map](figs/gibuu_eff_map_muon_p_costheta.png)
![muon_p_costheta projections](figs/gibuu_eff_proj_muon_p_costheta.png)

### ε(leading proton p, cos theta_p)

| leading proton p \ cos theta_p | [0.342, 0.4243) | [0.4243, 0.5065) | [0.5065, 0.5888) | [0.5888, 0.671) | [0.671, 0.7533) | [0.7533, 0.8355) | [0.8355, 0.9178) | [0.9178, 1) |
|---|---|---|---|---|---|---|---|---|
| [0.5, 0.625) | 0.114 ± 0.001 | 0.217 ± 0.001 | 0.296 ± 0.001 | 0.341 ± 0.002 | 0.350 ± 0.002 | 0.370 ± 0.002 | 0.377 ± 0.002 | 0.382 ± 0.002 |
| [0.625, 0.75) | 0.307 ± 0.001 | 0.377 ± 0.001 | 0.413 ± 0.002 | 0.432 ± 0.002 | 0.441 ± 0.002 | 0.452 ± 0.002 | 0.458 ± 0.002 | 0.428 ± 0.002 |
| [0.75, 0.875) | 0.283 ± 0.002 | 0.352 ± 0.001 | 0.386 ± 0.001 | 0.404 ± 0.002 | 0.408 ± 0.002 | 0.425 ± 0.002 | 0.427 ± 0.002 | 0.405 ± 0.002 |
| [0.875, 1) | 0.219 ± 0.002 | 0.255 ± 0.002 | 0.278 ± 0.001 | 0.290 ± 0.002 | 0.309 ± 0.002 | 0.332 ± 0.002 | 0.351 ± 0.002 | 0.351 ± 0.002 |
| [1, 1.1) | 0.153 ± 0.003 | 0.165 ± 0.002 | 0.178 ± 0.002 | 0.195 ± 0.002 | 0.215 ± 0.002 | 0.246 ± 0.003 | 0.274 ± 0.003 | 0.282 ± 0.003 |

![proton_p_costheta map](figs/gibuu_eff_map_proton_p_costheta.png)
![proton_p_costheta projections](figs/gibuu_eff_proj_proton_p_costheta.png)

| grid | ⟨ε⟩ | background at data POT |
|---|---|---|
| muon p | 0.3259 | 202110 |
| muon θ | 0.3259 | 202110 |
| muon p_T | 0.3259 | 202105 |
| leading proton p | 0.3259 | 170230 |
| leading proton θ | 0.3259 | 198410 |
| leading proton p_T | 0.3259 | 201926 |
| δp_T | 0.3259 | 201827 |
| δp_T (fine) | 0.3259 | 201827 |
| δp_Tx | 0.3259 | 201982 |
| δp_Ty | 0.3259 | 201982 |
| δα_T | 0.3259 | 201982 |
| φ_T | 0.3259 | 201982 |
| δp_L | 0.3259 | 201947 |
| p_n | 0.3259 | 201982 |
| cos θ_μ | 0.3259 | 202110 |
| cos θ_p | 0.3259 | 198410 |
| muon p × cos θ_μ | 0.3259 | 202110 |
| proton p × cos θ_p | 0.3259 | 166633 |

## 3. The factorised ansatz w = ε_μ(p_μ, cos θ_μ) · ε_p(p_p, cos θ_p) / ⟨ε⟩ — closure on the MC

Weighting the MC's own truth signal with the maps and comparing with the actually selected signal per true bin: total ratio 1.0032. Per grid (the error the shortcut makes on that grid):

| grid | ansatz / actual | max |rel. dev.| (bins > 50 events) |
|---|---|---|
| muon p | 1.0032 | 0.043 |
| muon θ | 1.0032 | 0.054 |
| muon p_T | 1.0032 | 0.237 |
| leading proton p | 1.0032 | 0.088 |
| leading proton θ | 1.0032 | 0.070 |
| leading proton p_T | 1.0032 | 0.061 |
| δp_T | 1.0032 | 0.458 |
| δp_T (fine) | 1.0032 | 0.398 |
| δp_Tx | 1.0032 | 0.124 |
| δp_Ty | 1.0032 | 0.507 |
| δα_T | 1.0032 | 0.039 |
| φ_T | 1.0032 | 0.164 |
| δp_L | 1.0032 | 0.458 |
| p_n | 1.0032 | 0.582 |
| cos θ_μ | 1.0032 | 0.022 |
| cos θ_p | 1.0032 | 0.012 |
| muon p × cos θ_μ | 1.0032 | 0.098 |
| proton p × cos θ_p | 1.0032 | 0.133 |

![ansatz dpt](figs/gibuu_ansatz_dpt.png)
![ansatz alpha](figs/gibuu_ansatz_alpha.png)
![ansatz pn](figs/gibuu_ansatz_pn.png)
![ansatz muon_theta](figs/gibuu_ansatz_muon_theta.png)
![ansatz proton_theta](figs/gibuu_ansatz_proton_theta.png)

## 4. Background of the selected sample (official MC, at the data POT)

Selected reco events that are not truth signal inside the fiducial volume, by the category of the reco rows' truth, plus feed-in (signal whose true observable lies outside the grid). The MC is the unweighted central value.

| grid | total | single π± | single π⁰ | multi-π | no pion | feed-in |
|---|---|---|---|---|---|---|
| muon p | 202110 | 94927 | 45622 | 25375 | 36185 | 0 |
| muon θ | 202110 | 94927 | 45622 | 25375 | 36185 | 0 |
| muon p_T | 202105 | 94925 | 45622 | 25374 | 36184 | 0 |
| leading proton p | 170230 | 84202 | 40297 | 22200 | 23531 | 0 |
| leading proton θ | 198410 | 94477 | 45456 | 25320 | 33157 | 0 |
| leading proton p_T | 201926 | 94870 | 45603 | 25367 | 36085 | 0 |
| δp_T | 201827 | 94840 | 45570 | 25334 | 36082 | 0 |
| δp_T (fine) | 201827 | 94840 | 45570 | 25334 | 36082 | 0 |
| δp_Tx | 201982 | 94892 | 45609 | 25371 | 36111 | 0 |
| δp_Ty | 201982 | 94892 | 45609 | 25371 | 36111 | 0 |
| δα_T | 201982 | 94892 | 45609 | 25371 | 36111 | 0 |
| φ_T | 201982 | 94892 | 45609 | 25371 | 36111 | 0 |
| δp_L | 201947 | 94879 | 45600 | 25366 | 36103 | 0 |
| p_n | 201982 | 94891 | 45609 | 25371 | 36110 | 0 |
| cos θ_μ | 202110 | 94927 | 45622 | 25375 | 36185 | 0 |
| cos θ_p | 198410 | 94477 | 45456 | 25320 | 33157 | 0 |
| muon p × cos θ_μ | 202110 | 94927 | 45622 | 25375 | 36185 | 0 |
| proton p × cos θ_p | 166633 | 83791 | 40141 | 22151 | 20551 | 0 |

## 5. Data versus GiBUU signal + MC background

Three predictions per grid: **full** = GiBUU truth cells × efficiency × migration (the binned response learned on the official MC) + background; **eff-only** = GiBUU truth cells × efficiency, no migration; **ansatz** = GiBUU truth events × w from the two maps, histogrammed in truth bins. −2lnL is the Baker–Cousins Poisson likelihood ratio of the data against the full prediction.

| grid | data | full | eff-only | ansatz | background | data/full | −2lnL/ndf (full) |
|---|---|---|---|---|---|---|---|
| muon p | 329653 | 377228 | 377228 | 380660 | 202110 | 0.874 | 7676.8/8 |
| muon θ | 329653 | 376281 | 376281 | 380660 | 202110 | 0.876 | 7362.9/14 |
| muon p_T | 329650 | 373754 | 373759 | 380655 | 202105 | 0.882 | 8770.4/10 |
| leading proton p | 289289 | 334950 | 343748 | 348781 | 170230 | 0.864 | 7240.4/5 |
| leading proton θ | 323116 | 374157 | 375129 | 376960 | 198410 | 0.864 | 12692.1/10 |
| leading proton p_T | 329382 | 377238 | 377278 | 380476 | 201926 | 0.873 | 7593.0/7 |
| δp_T | 329293 | 377335 | 377429 | 380377 | 201827 | 0.873 | 8589.6/6 |
| δp_T (fine) | 329293 | 377708 | 377801 | 380377 | 201827 | 0.872 | 9125.5/9 |
| δp_Tx | 329467 | 376210 | 376237 | 380533 | 201982 | 0.876 | 6932.7/8 |
| δp_Ty | 329467 | 376260 | 376287 | 380533 | 201982 | 0.876 | 6895.0/11 |
| δα_T | 329467 | 375691 | 375718 | 380533 | 201982 | 0.877 | 6144.1/5 |
| φ_T | 329467 | 376985 | 377012 | 380533 | 201982 | 0.874 | 8067.1/6 |
| δp_L | 329420 | 379475 | 379521 | 380498 | 201947 | 0.868 | 14462.6/7 |
| p_n | 329467 | 379275 | 379302 | 380532 | 201982 | 0.869 | 13069.6/10 |
| cos θ_μ | 329653 | 376542 | 376542 | 380660 | 202110 | 0.875 | 7047.5/8 |
| cos θ_p | 323116 | 373895 | 374886 | 376960 | 198410 | 0.864 | 12829.8/8 |
| muon p × cos θ_μ | 329653 | 376299 | 376299 | 380660 | 202110 | 0.876 | 9875.0/64 |
| proton p × cos θ_p | 282881 | 333725 | 344013 | 345184 | 166633 | 0.848 | 15331.9/40 |

![muon_p](figs/gibuu_overlay_muon_p.png)
*Figure 1 — muon p.*

![muon_theta](figs/gibuu_overlay_muon_theta.png)
*Figure 2 — muon θ.*

![muon_pt](figs/gibuu_overlay_muon_pt.png)
*Figure 3 — muon p_T.*

![proton_p](figs/gibuu_overlay_proton_p.png)
*Figure 4 — leading proton p.*

![proton_theta](figs/gibuu_overlay_proton_theta.png)
*Figure 5 — leading proton θ.*

![proton_pt](figs/gibuu_overlay_proton_pt.png)
*Figure 6 — leading proton p_T.*

![dpt](figs/gibuu_overlay_dpt.png)
*Figure 7 — δp_T.*

![dpt_fine](figs/gibuu_overlay_dpt_fine.png)
*Figure 8 — δp_T (fine).*

![dptx](figs/gibuu_overlay_dptx.png)
*Figure 9 — δp_Tx.*

![dpty](figs/gibuu_overlay_dpty.png)
*Figure 10 — δp_Ty.*

![alpha](figs/gibuu_overlay_alpha.png)
*Figure 11 — δα_T.*

![phi](figs/gibuu_overlay_phi.png)
*Figure 12 — φ_T.*

![pl](figs/gibuu_overlay_pl.png)
*Figure 13 — δp_L.*

![pn](figs/gibuu_overlay_pn.png)
*Figure 14 — p_n.*

![muon_costheta](figs/gibuu_overlay_muon_costheta.png)
*Figure 15 — cos θ_μ.*

![proton_costheta](figs/gibuu_overlay_proton_costheta.png)
*Figure 16 — cos θ_p.*

![muon_p_costheta](figs/gibuu_folded_muon_p_costheta.png)
*Figure 17 — muon p × cos θ_μ.*

![proton_p_costheta](figs/gibuu_folded_proton_p_costheta.png)
*Figure 18 — proton p × cos θ_p.*

## 6. Caveats

- The background and the efficiency come from the unweighted central-value official MC (no MINERvA tune, flux or detector weights); the data/MC ratio of the selected sample is 0.84 with that MC, so absolute agreement or disagreement of the overlay carries that uncertainty.
- The fiducial volume and n_nucleons are the inclusive channel's (the paper's own are unstated); the flux integral is the inclusive channel's Φ.
- GiBUU: carbon only (hydrogen gives no signal), the channel flux rebinned to 0.5 GeV, card physics from GiBUU's MINERvA-ME card, statistics and weights as recorded in the sample's `gibuu_run.json`.
- The ansatz prediction has no detector migration; its closure table (Sec. 3) is the size of that approximation on each grid.

