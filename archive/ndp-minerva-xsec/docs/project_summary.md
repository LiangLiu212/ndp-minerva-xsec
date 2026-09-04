# Project summary — reproducing arXiv:2106.16210

Physics summary of this project: the MINERvA medium-energy (ME) FHC **inclusive
charged-current ν_μ** double-differential cross section **d²σ/dp_T dp_∥** on
hydrocarbon, ⟨E_ν⟩ ≈ 6 GeV (Ruterbories *et al.*, Phys. Rev. D **104**, 092007),
reproduced from the MINERvA Open Data MasterAnaDev AnaTuples in pure Python.

Pipeline: `ndp-minerva-xsec`, four stages (inputs → cuts/signal-def →
ingredients → unfolding+normalization; see `README.md`). **Current scope: playlist
1A only, stat-only central value; the 12-playlist combine is the open gate.**

**Headline:** the full chain reproduces the published cross section to **0.13%
integrated (205/205 reported cells)**; the assembled systematic budget is at
**88.5%** of the paper's. Paper analysis-only extract:
`../../ndp-minerva-data-release-exploration/papers/minerva/2106.16210/paper_2106.16210.md`.

---

## 1. Cross-section formula

The reported double-differential cross section, per (p_T, p_∥) cell *(i, j)*:

$$\frac{d^{2}\sigma}{dp_{T}\,dp_{\parallel}}\bigg|_{ij} = \frac{1}{\Phi \cdot N \cdot \Delta p_{T,i}\,\Delta p_{\parallel,j}} \cdot \frac{U_{ij}\big[N^{data}_{ij}-N^{bkg}_{ij}\big]}{\varepsilon_{ij}}$$

| Term | Meaning | Paper | Ours |
|---|---|---|---|
| N_data − N_bkg | background-subtracted reco counts (bkg subtracted bin-by-bin) | bkg 0.2% | — |
| U | D'Agostini unfolding operator (reco → true) | 10 iterations | self-closure 6e-11 |
| ε | MC efficiency-acceptance correction | 64% | 0.657 |
| Φ | integrated flux × exposure | 6.32×10⁻⁸ cm⁻² POT⁻¹ × 10.61×10²⁰ POT | Φ_int 6.23×10⁻⁸ |
| N | target nucleons (fiducial) | 3.23×10³⁰ ± 1.4% | 3.2353×10³⁰ (+0.16%) |
| Δp_T, Δp_∥ | bin widths | — | — |

1D projections are derived **from the 2D result** (integrate the other axis within
its reported range), not measured independently. **Our integrated σ = 3.035×10⁻³⁸
cm²/nucleon, ratio to published 0.9987.** Implementation: `extract_xsec.py`
(`xsec/flux.py`, `xsec/targets.py`, `xsec/binning.py`).

---

## 2. The input

- **Beam / detector (paper):** NuMI FHC, ⟨E_ν⟩ ≈ 6 GeV; MINERvA fine-grained
  scintillator tracker (fiducial **5.48 t** hydrocarbon) + MINOS near detector for
  muon charge sign and momentum; exposure **10.61×10²⁰ POT** (Sep 2013 – Feb 2017).
- **Simulation:** GENIE 2.12.6 + **MINERvA Tune v1** (Valencia RPA on QE, 2p2h with
  +50% empirical enhancement, non-resonant π reduced 43%); PPFX NuMI flux with the
  in-situ ν-e scattering constraint.
- **Our data source:** MINERvA Open Data, ME FHC playlists **1A–1P**, MasterAnaDev
  AnaTuples (Data + StandardMC), **streamed via xrootd** with branch-pruned `uproot`
  — the heavy tuples are never downloaded; dataset specs + streamed fingerprints in
  `config/datasets/*.json` (`summarize_inputs.py`). POT ledger Σ = **10.574×10²⁰ =
  99.66%** of the paper's exposure (open data carries the full exposure).
- The one local-fetch exception is the flux/reweight-files tarball (CV weight stage).
- **Reference implementation:** the MINERvA-101 cross-section tutorial (same tuples,
  CC-inclusive chain), extended here to 2D and the full ME phase space.

---

## 3. The cuts

Reconstruction-level selection, applied identically to data and MC — deliberately
few, because the measurement is fully inclusive (`xsec/cuts.py`, `xsec/kinematics.py`):

1. Muon track **matched MINERvA ↔ MINOS**.
2. Muon charge **negative** (μ⁻, from MINOS magnetization).
3. Reconstructed vertex inside the **fiducial volume** (850 mm apothem hexagon, 5.48 t).
4. **θ_μ < 20°**, **1.5 ≤ p_∥ ≤ 60 GeV/c**, **p_T < 4.5 GeV/c** — mirrors the
   reporting phase space and enforces MINOS acceptance (note the induced geometric
   limit Q² ≲ p_∥²/8).

**Selected data:** 4,105,696 events. **Predicted background:** 8,655 (0.2%),
subtracted bin-by-bin before unfolding (low-p_∥/low-p_T = NC π→μ fakes; high-p_T =
ν̄_μ charge-misID). **Average efficiency 64%**, driven by MINERvA–MINOS geometric
acceptance, not hadronic cuts. The **signal definition** (`xsec/signal.py`) is the
same phase space at truth level. Our selector is certified by cut-flow parity
(844 / 43643 / 43539 / 104) against the exploration repo's golden selector.

---

## 4. The unfolding

- **Method:** D'Agostini iterative unfolding (RooUnfold in the paper; `xsec/unfold.py`
  here), **10 iterations** — fixed via a χ²-vs-truth scan over 10 pseudodata
  variations (2p2h strength, QE RPA, non-res π, res π reweights, plus a data-driven
  reweight); the required iterations never exceeded 10.
- **Binning:** **16 p_∥ × 14 p_T = 224 cells** (205 reported; 19 empty by
  acceptance). p_∥ edges 1.5 → 60 GeV/c; p_T edges 0 → 4.5 GeV/c (`xsec/binning.py`).
- **Migration:** reco ↔ true in (p_T, p_∥); we use count-conserving flat slots
  (288 = 16 × 18, under/overflow retained), oriented `[reco, true]` so the unfolder's
  efficiency numerator is free from the column sum.
- **Order:** background-subtract → unfold (U) → divide by efficiency (ε) → normalize
  by Φ · N · POT · bin width. MC self-closure 6e-11 for the unfolder; 2.9e-15 through
  the whole extraction chain.

---

## 5. Systematic uncertainties

The paper uses **three categories** — Flux / Detector / Interaction model — refined
to seven curves in Fig. 8. Per-cell median fractional uncertainty (playlist 1A):

| Fig-8 curve | Category | Our median | Status / validation |
|---|---|---:|---|
| **Muon reconstruction** | Detector | **3.79%** | done; energy-scale **0.99** vs `cov_energyscale` |
| **Flux** | Flux | **3.56%** | shape-resolved (100 PPFX); **0.88** + off-diag corr **0.85** vs `cov_flux` |
| Statistical (1A) | — | 4.70% | data-Poisson toys; **0.90** vs `cov_stat` — *inflated, 1A only* |
| Models | Interaction | ~1.3% | GENIE 1.19 + 2p2h 0.31 + RPA 0.10 + GenieRvx1pi 0.03 |
| Normalization | Flux | 1.40% | exact flat band (target-N count); trivial rank-1 |
| Hadronic Response | Detector | — | **not built** (~1%) |

- **Construction:** each group is built from per-universe re-extractions of the full
  cross section, summed as covariances (`xsec/systematics.py`, `make_*_universes.py`,
  `assemble_*.py`); muon-reco is one reco-only pass over 8 bands / 16 universes.
- **Totals:** systematic-only **5.54% = 88.5%** of the published 6.26%; total **7.51%**
  vs published 6.83%. The total overshoots **only** because the 1A statistical band
  (4.7%) stands in for the full-dataset stat (~1.5%) — the fair comparison is the
  systematic-only number.
- **Anc-validated** on all three dedicated published covariance files: flux shape,
  statistical, muon energy scale.
- **Fig-8-format breakdown:** `results/img/total_1A/fig8_errors_{pt,pz}.png`
  (`fig8_uncertainties.py`).
- **Remaining:** Hadronic Response (geant4/response), the "new" GENIE bands
  (MaRES⊗NormCCRES + FaCCQE), the (trivial) normalization band, the MINOS-band
  flux-weight correlation, and the **12-playlist combine** — the gate for a clean
  total-covariance validation.

---

## Pointers

- Code workflow + stage scripts: `../README.md`.
- Central-value result + full budget: `results/2026-06-16_xsec_1A_with_systematics.md`.
- Systematics status: `results/2026-06-16_systematics_summary.md`.
- Per-stage detail: `results/2026-06-16_{flux_shape,muon_reconstruction,stat,
  energyscale,2p2h,rpa,geniervx1pi,total}_*.md`.
- Paper analysis-only extract:
  `../../ndp-minerva-data-release-exploration/papers/minerva/2106.16210/paper_2106.16210.md`.
