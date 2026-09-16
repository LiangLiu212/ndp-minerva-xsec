# GiBUU against the MINERvA 1μ1p data in reconstruction space

**What this is.** A GiBUU 2025 prediction for the CCQE-like νμ + leading-proton sample of
arXiv:2503.15047, pushed forward into the reconstructed quantities the experiment actually measured
and compared bin by bin with the full medium-energy FHC Open Data. The data are never unfolded. The
generator's truth events are normalised absolutely (cross section × flux × POT × nucleons), folded
through the detector response learned from the official MINERvA MC, and added to the background
predicted by that same MC. The comparison covers all 18 grids: muon p, θ, cos θ, p_T; leading-proton
p, θ, cos θ, p_T; and the transverse-imbalance set δp_T (two binnings), δp_Tx, δp_Ty, δα_T, φ_T,
δp_L, p_n; plus the two two-dimensional efficiency grids.

**Status:** 2026-09-16. First GiBUU comparison on this channel. The background and the detector
response come from the unweighted central-value official MC, with no MINERvA tune, flux constraint or
systematic universes, so the absolute agreement carries that uncertainty; shapes are the robust part.

**Provenance.** Sample `runs/_generator_cache/gibuu_1f8bcb950be1dcd0/` (campaign
`grid/campaigns/gibuu_me_c12_scan_2026-09/`), comparison run
`runs/2026-09-16_gibuu_2025_me_fhc_c12__minerva_me_ccqelike_1mu1p__all_3/`, efficiency and background
from `runs/2026-09-15_efficiency_minerva_me_ccqelike_1mu1p/`. Machine-rendered tables and all 18
overlays: `report/GiBUU_overlay_FHC.md`. Every number below is quoted from those directories.

---

## 1. Inputs

| item | value |
|---|---|
| data | ME FHC playlists 1A–1P, 1.0574 × 10²¹ POT_Used, 329 653 selected events |
| official MC | 4.9784 × 10²¹ POT; supplies the efficiency, the migration and the background |
| GiBUU sample | GiBUU 2025 patch 5, νμ CC on ¹²C, 979 395 events, 153 grid jobs |
| GiBUU σ_CC | 4.1245 × 10⁻³⁸ cm² per nucleon (flux-averaged over the covered flux) |
| normalisation | N = σ_cell × 3.23 × 10³⁰ nucleons × 6.32 × 10⁻⁸ cm⁻² POT⁻¹ × 1.0574 × 10²¹ POT |
| selection | `minerva_ccqelike_1mu1p_v0`, identical for data, MC and the folded prediction |

## 2. How the GiBUU sample was made, and why that way

The obvious route fails. GiBUU's flux-averaged event mode (`nuXsectionMode = 16`) samples the lepton
kinematics uniformly and compensates with event weights. Under the NuMI medium-energy flux those
weights are extremely uneven for exactly the channel this analysis selects: in a 1000-ensemble pilot,
1 140 quasi-elastic events carried the statistical power of 30, and a 328-event 1μ1p signal sample
the power of 21. Over half of the quasi-elastic cross section sat in events more than twenty times
heavier than the mean. GiBUU's rejection mode, which would flatten the weights, aborts the run on any
weight above its ceiling, and the quasi-elastic tail exceeded ceilings of 500 and 1 500 within a few
hundred test nucleons.

The sample was therefore generated as a **flux-weighted energy scan**: one fixed-energy
`integratedSigma` run per 0.5 GeV point from 2 to 60 GeV, 116 points in 156 grid jobs. In that mode
GiBUU integrates the cross section for each test nucleon and draws the kinematics by rejection
internally, so the event weights are uniform to about 35 % and one event is produced per test
nucleon. The points are recombined with each point's share of the true flux,

> weight = perweight × (flux fraction of the point's bin) / (number of jobs at that point),

so the merged weights sum to the flux-averaged cross section per nucleon over the covered range and
the per-point σ_CC(E) is available as a by-product. Ensembles were allocated in proportion to flux × E
so the merged sample is not dominated by a handful of energies.

**Validation.** Each job's summed weights must equal the cross section GiBUU writes to its own
absorption file; the merge refuses a job that disagrees by more than one part in a thousand. The
per-point σ_CC(E) tracks the GENIE G18_02a spline for the same target closely up to 30 GeV and falls
a few per cent below it above. Two independent seeds at 6 GeV agreed on σ_QE to 1 %.

![sigma of energy](figs/gibuu_sigma_of_energy.png)
*Figure 1 — GiBUU's cross section per energy point against the GENIE spline, with each point's flux
weight and event count.*

**Coverage.** 115 of the 116 points are present, carrying 94.6 % of the 0–100 GeV flux. The one
missing point is 2.25 GeV, where the signal requires a muon above 2 GeV/c together with a proton of
0.5 to 1.1 GeV/c, which is kinematically closed. The prediction is therefore complete for the signal
even though the flux bookkeeping is not.

## 3. From truth to reconstructed events

Three steps, all taken from the official MC and applied identically to any model:

1. **Absolute normalisation.** GiBUU's per-nucleon cross section per truth cell is multiplied by the
   nucleon count of the fiducial volume, the flux integral per POT and the data POT. The factor is
   2.158 × 10⁶ events per 10⁻³⁸ cm².
2. **Efficiency and migration.** The truth cells are multiplied by the selection efficiency of that
   cell and redistributed over reconstructed cells by the migration matrix, both learned from the
   paired official MC on each grid.
3. **Background.** The selected MC events that are not truth signal in the fiducial volume are added,
   scaled by the POT ratio: 201 827 events, composed of 94 927 single-π±, 45 622 single-π⁰, 25 375
   multi-pion and 36 185 no-pion events.

The report also carries two cruder variants for reference: **efficiency only**, which skips the
migration, and the **factorised ansatz**, which weights each truth event by ε_μ(p_μ, cos θ_μ) ·
ε_p(p_p, cos θ_p) / ⟨ε⟩ from the two maps. They differ from the full fold by well under a per cent on
the muon variables and by up to 3 % on the proton momentum, which is the size of the migration
correction there.

**A consistency check of the whole normalisation chain.** GiBUU predicts 533 063 truth signal events
in the fiducial volume at the data exposure, computed from a cross section, a flux integral and a
nucleon count. The official MC, normalised by POT alone and never touching those three numbers,
contains 580 715. Two entirely different normalisation routes agree to 8 %, which is within the flux
and nucleon-count uncertainties and is the strongest evidence that the absolute chain is right.

## 4. What the comparison shows

### 4.1 Rates

| | data | signal predicted | background | total | data / prediction |
|---|---|---|---|---|---|
| GiBUU + MC background | 329 293 | 175 508 | 201 827 | 377 335 | **0.873** |
| GENIE central value | 329 653 | 189 275 | 202 110 | 391 384 | **0.842** |

Both generators predict more selected events than the data. GiBUU predicts 7 % fewer signal events
than GENIE and lands about 4 % closer to the data. Across the 18 grids GiBUU's ratio sits between
0.848 and 0.882, the spread reflecting the different acceptance of each grid rather than a physics
difference.

Neither number should be read as a measurement of a cross section. The background and the response
come from an MC with no tune and no flux constraint, and that MC is itself 16 % above the data. The
honest statement is that GiBUU and GENIE bracket the same 12 to 16 % overshoot, with GiBUU slightly
lower.

### 4.2 Shapes, which is where the two models genuinely differ

Removing the normalisation from each generator separately and looking at what is left:

| grid | GiBUU shape rms | GiBUU worst bin | GENIE shape rms | GENIE worst bin |
|---|---|---|---|---|
| δα_T | **0.022** | 0.043 | 0.059 | 0.135 |
| leading proton p | 0.051 | 0.083 | 0.054 | 0.094 |
| δp_Tx | **0.061** | 0.151 | 0.097 | 0.271 |
| δp_T | 0.072 | 0.145 | **0.018** | 0.042 |
| φ_T | 0.071 | 0.119 | 0.065 | 0.115 |
| δp_Ty | 0.101 | 0.241 | 0.120 | 0.335 |
| muon θ | **0.105** | 0.315 | 0.132 | 0.401 |
| proton p_T | 0.117 | 0.291 | 0.111 | 0.317 |
| proton θ | 0.131 | 0.332 | **0.064** | 0.171 |
| muon p | 0.148 | 0.386 | 0.138 | 0.316 |
| muon p_T | **0.165** | 0.416 | 0.203 | 0.485 |
| δp_L | 0.194 | 0.341 | 0.151 | 0.252 |
| p_n | 0.194 | 0.519 | **0.070** | 0.133 |

GiBUU describes the angular correlation δα_T markedly better than GENIE (2 % against 6 %) and is also
better on δp_Tx, muon θ and muon p_T. GENIE describes δp_T, p_n and the proton angle better.

**The low-Q² region.** Both models over-predict forward, low-transverse-momentum muons, the classic
signature of the missing low-Q² suppression and 2p2h enhancement that the MINERvA tune supplies.
GiBUU's deficit is less severe: data over prediction in the first muon p_T bin is 0.52 against GENIE's
0.43, and in the first muon θ bin 0.60 against 0.50. Neither model is close there without the tune.

**The Fermi-motion peak.** This is the sharpest difference. At p_n below 0.1 GeV/c the data stand 32 %
*above* the GiBUU prediction, while they sit 27 % *below* GENIE's. The same reversal appears in the
first δp_T bin, where GiBUU's ratio is 1.00 and GENIE's 0.85, and in δp_L below 0.1 GeV/c. GiBUU's
nuclear ground state and its treatment of final-state interactions move strength out of the
quasi-elastic peak that GENIE leaves in, and the data prefer something in between, closer to GENIE at
the peak and closer to GiBUU in the tails.

![p_n overlay](figs/gibuu_overlay_pn.png)
*Figure 2 — p_n: data, GiBUU signal folded onto the background, and the two cruder prediction
variants. The ratio panel shows the reversal at low p_n.*

![δp_T overlay](figs/gibuu_overlay_dpt.png)
*Figure 3 — δp_T. GiBUU reproduces the first bin exactly and under-predicts the 0.2 to 0.6 GeV/c
region by 18 %.*

![δα_T overlay](figs/gibuu_overlay_alpha.png)
*Figure 4 — δα_T, the variable GiBUU describes best: a flat 13 % normalisation offset and no shape
trend.*

![muon θ overlay](figs/gibuu_overlay_muon_theta.png)
*Figure 5 — muon angle, where both generators show the low-Q² trend and GiBUU is the milder of the two.*

![proton momentum overlay](figs/gibuu_overlay_proton_p.png)
*Figure 6 — leading-proton momentum, the grid where the migration correction matters most (the dashed
efficiency-only curve differs from the full fold by 3 %).*

### 4.3 Where the difference comes from

The two generators build the same signal out of different processes:

| process | GiBUU share of the signal cross section | GENIE share of the fiducial signal |
|---|---|---|
| quasi-elastic | 36.7 % | 46.7 % |
| resonance with the pion absorbed | 31.2 % | 25.5 % |
| 2p2h | 16.5 % | 21.0 % |
| deep inelastic | 15.6 % | 6.9 % |

GiBUU reaches a similar total with a third less quasi-elastic, more resonance production and more
than twice the deep-inelastic contribution. Since the quasi-elastic component is what populates the
low-p_n, low-δp_T peak, that redistribution is exactly what produces the shape differences of
Section 4.2. GiBUU's signal cross section is 0.2470 × 10⁻³⁸ cm² per nucleon, 6.0 % of its charged-current
total.

### 4.4 Goodness of fit

The report quotes Baker–Cousins −2lnλ per grid. With 330 000 selected events the statistical errors
are a few tenths of a per cent, so a 13 % normalisation offset dominates every number: values run
from 6 100 for δα_T over five bins to 15 300 for the two-dimensional proton grid over forty. These
numbers are not usable as a fit quality until the MC tune and flux weights are applied and their
systematic covariance is included. The shape table of Section 4.2 is the meaningful comparison today.

## 5. Caveats

- **The background and the response are unweighted central-value MC.** No MINERvA tune, no flux
  constraint, no detector systematics. That MC over-predicts the selected data by 16 % on its own, so
  the 13 % offset seen here is not a GiBUU result alone.
- **The efficiency has a playlist dependence** not yet understood: it falls from 0.35 in playlists
  1A to 1E to 0.28 in 1O and 1P at constant purity. The average is used throughout.
- **The fiducial volume and nucleon count** are the inclusive channel's, because the paper does not
  state its own. The flux integral is likewise inherited.
- **GiBUU is carbon only.** Hydrogen, 8.2 % of the CH target by mass, cannot produce a 1μ1p signal
  event, so the carbon sample scaled by all nucleons is the right per-nucleon treatment, but the
  small nuclear-target differences of the real tracker are not modelled.
- **Discretising the flux in 0.5 GeV steps** is a modelling choice of this generation method. Inside
  a step the neutrino energy is fixed, so the muon spectrum is slightly staircased at a scale far
  below any bin used here.
- **The 2.25 GeV point is missing** from the flux bookkeeping, 2.7 % of the flux, with no signal.
- **The purity of the selection is 0.48 against the paper's 0.60**, an open question of the
  selection transcription; a different background level would move the ratios in Section 4.1.

## 6. Reproduce

```bash
cd /exp/dune/data/users/liangliu/ndp-dev/ndp-platform
# the sample (grid/README.md has the full recipe)
python -m ndp gibuu plan gibuu_me_c12_scan_2026-09 models/gibuu_2025_me_fhc_c12.yaml --channel minerva_me_ccqelike_1mu1p
grid/stage_gibuu_payload.sh runs/_generator_cache/gibuu_<fingerprint>
python -m ndp gibuu submit-cmd gibuu_me_c12_scan_2026-09 --tar-label ndp-gibuu-v1   # then status / harvest
python -m ndp gibuu merge gibuu_me_c12_scan_2026-09 --channel minerva_me_ccqelike_1mu1p
# the response, efficiency and background
python -m ndp surrogate build --channel minerva_me_ccqelike_1mu1p --measurement all --kind binned
python -m ndp efficiency run --channel minerva_me_ccqelike_1mu1p
# the comparison and the reports
python -m ndp run models/gibuu_2025_me_fhc_c12.yaml --channel minerva_me_ccqelike_1mu1p \
    --measurement all --modes folded --efficiency-run runs/2026-09-15_efficiency_minerva_me_ccqelike_1mu1p
python report/make_overlay_report.py --efficiency runs/2026-09-15_efficiency_minerva_me_ccqelike_1mu1p \
    --overlay runs/<the run above> --out report/GiBUU_overlay_FHC.md
```

## 7. What would sharpen this

1. Apply the MINERvA tune and flux weights to the MC, so that the background and the response stop
   carrying a 16 % normalisation error into the comparison.
2. Settle the purity gap against the paper, which sets the background level.
3. Run the same energy scan for the alternative GiBUU configurations, in particular with the
   in-medium widths off, to see how much of the δα_T and p_n behaviour is transport and how much is
   the initial state.
4. Add the systematic universes to the folded prediction so the goodness-of-fit numbers mean
   something.
