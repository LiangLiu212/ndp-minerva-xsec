# MINERvA CCQE-like 1μ1p / TKI analysis on the open data — playlist 1A

**Status:** rendered 2026-09-15 from runs of 2026-09-14. Truth-level signal decided (paper definition); reconstruction-level selection transcribed from the paper with the unpublished cut values marked as defaults; data versus the official MC at reconstruction level, shape only (the MC is the unweighted central value).

**Reference analysis:** J. Kleykamp et al. (MINERvA), *Measurement of the A dependence of the νμ charged-current quasielastic-like cross section as a function of muon and proton kinematics at ⟨Eν⟩ ∼ 6 GeV*, arXiv:2503.15047, Phys. Rev. D 112, 052005 (2025). Only the CH (tracker) result is reproduced here.

**Code:** `ndp-platform` (report rendered at commit `cbad5f7`; signal run at `50f1099`, selection run at `75728d6`). Channel manifest `channels/minerva_me_ccqelike_1mu1p.yaml`; signal `ndp/channels/signal.py`; observables `ndp/channels/observables.py`, `ndp/channels/reco_observables.py`; selection `ndp/adapters/minerva_anatuple.py::ccqelike_1mu1p_cutflow` via `ndp/channels/selections.py`; diagnostics `ndp/diagnostics.py` (`python -m ndp signal`, `python -m ndp selection`); grid processing `ndp/grid/`, playlist products `ndp/products.py`.

**Runs (every number below is quoted from these directories):** `runs/2026-09-14_signal_minerva_me_ccqelike_1mu1p_FHC1A/` (truth level) and `runs/2026-09-14_selection_minerva_me_ccqelike_1mu1p_FHC1A/` (reconstruction level). This file was rendered by `report/make_report.py` from them.

---

## 1. Inputs

| item | value | source |
|---|---|---|
| Data | ME FHC playlist 1A: 253 MasterAnaDev AnaTuple files, streamed on the FNAL grid (campaign `grid/campaigns/fhc_2026-09`) | MINERvA Open Data |
| Data exposure | 8.969 × 10¹⁹ POT | Meta tree `POT_Used`, summed over files (`pot_<pl>_data.json`) |
| Official MC | StandardMC, same playlist: 41 files, GENIE 2.12.6, unweighted central value (no MINERvA-tune, flux or GENIE weights) | MINERvA Open Data |
| MC exposure | 4.072 × 10²⁰ POT; data/MC POT scale 0.22028 | Meta tree, summed |
| Truth (MC) | 5 083 107 CC events in the skim box (z 5880–8522 mm, apothem 900 mm), 4 209 798 inside the tracker fiducial | `truth_<pl>_skim.npz` |
| Reco candidates (MC) | 7 605 892 reconstructed candidates, with the reco rows' truth | `reco_<pl>_mc.npz`, `reco_<pl>_truthcols_skim.npz` (cache version 3) |
| Reco candidates (data) | 2 791 649 reconstructed candidates | `reco_<pl>_data.npz` |
| Paper exposure | 10.61 × 10²⁰ POT, 218 000 selected CH events | arXiv:2503.15047 |

This data set is 8.5 % of the paper's exposure. The statistical error on the selected data/MC ratio is ±0.005; every ratio below is a shape statement because the MC carries no weights.

## 2. Signal definition (truth level)

### 2.1 The paper's definition

Verbatim from the paper's tex (`NuclTargTKI_PRD.tex` lines 293–301):

> "The interactions are considered signal if they have a muon with an angle with respect to the
> beam of <17° and a momentum within the range 2 GeV/c < p_μ < 20 GeV/c, and a proton with an
> angle <70° and a momentum in the range 500 MeV/c < p_p < 1100 MeV/c. Additionally, the
> interaction must not have mesons, baryons heavier than neutrons, or photons above 10 MeV.
> [...] For events with more than one proton matching the constraints, the highest momentum
> matching proton is used. [...] The photons with an energy less than 10 MeV are accepted since
> they can come from nuclear de-excitations."

As implemented (channel manifest, `signal:` block, status **decided**):

| item | requirement |
|---|---|
| interaction | νμ charged current (`nu_pdg = 14`, `current = CC`, primary lepton `lep_pdg = 13`) |
| muon | θμ < 17° with respect to the beam; 2 < pμ < 20 GeV/c |
| proton | ≥ 1 proton with θp < 70° and 0.5 < pp < 1.1 GeV/c; the highest-momentum proton inside this window is the **leading proton** used by every proton and TKI observable |
| vetoes | no mesons; no baryons heavier than the neutron; no photons above 10 MeV |
| allowed | any number of neutrons, photons ≤ 10 MeV, protons outside the window, nuclear remnants |
| phase space | true vertex inside the tracker fiducial box of the inclusive channel (z 5980–8422 mm, apothem 850 mm) — **open**: the paper does not state its CH fiducial volume |

Implementation choices, logged in `docs/decisions.md` (2026-09-14):

- Final-state classes are PDG **ranges**, not an enumerated list: meson = 100 ≤ |pdg| < 1000; heavy baryon = 1000 ≤ |pdg| < 10000 except p and n; antinucleons are *not* heavier than the neutron and are not vetoed; nuclear remnants (|pdg| ≥ 10⁹) and GENIE bookkeeping particles (2000000101) are ignored. Against MAT's enumerated `IsQELike` this differs on 6 of 397 604 νμ CC events of the certification file (2 anti-Λ, 1 anti-K⁰ events vetoed here but not by MAT; 3–4 events with a second final-state muon kept here, rejected by MAT).
- Both the muon and the proton angles are evaluated in the **NuMI beam frame** (rotation of the detector-frame momenta about x by −0.05887 rad, MAT's convention), the frame the paper measures angles in and the platform's default with the 2026-09-04 evidence.
- On the grid the definition is evaluated per file by the same code and stored as derived columns of the truth skim (`signal_minerva_ccqelike_1mu1p`, `lp_*`, veto counts); the numbers below re-evaluate it from those columns.

### 2.2 Truth-level result on the official MC

Cumulative cutflow on the truth skims (`cutflow.json` of the signal run):

| step | all skimmed events | inside tracker fiducial |
|---|---|---|
| all cached truth events | 5 083 107 | 4 209 798 |
| νμ CC with μ⁻ | 4 921 173 | 4 075 640 |
| muon window | 2 800 145 | 2 319 685 |
| ≥ 1 proton in window | 961 529 | 796 255 |
| no mesons | 278 050 | 230 244 |
| no heavy baryons | 278 047 | 230 241 |
| no photons > 10 MeV | **270 020** | **223 595** |

Signal fraction of fiducial νμ CC events: 5.49 % (223 595 of 4 075 640). Composition of the fiducial signal by GENIE process: QE 104 114 (46.6 %), RES with the pion absorbed 57 047 (25.5 %), 2p2h 46 875 (21.0 %), DIS 15 559 (7.0 %). Protons inside the window per signal event: one in 197 824 events, two in 22 614, three or more in 3 157.

Kinematics of the fiducial signal (median, 16th–84th percentile): muon p 5.08 (3.39–6.96) GeV/c, muon θ 7.2° (4.0–11.7°), leading proton p 0.741 (0.576–0.954) GeV/c, leading proton θ 53.5° (31.3–64.4°), leading proton pT 0.556 (0.363–0.752) GeV/c.

![truth-level signal kinematics](figs/pl1A_signal_kinematics.png)
*Figure 1 — Muon and leading-proton kinematics of the fiducial truth-level signal (beam frame), and the proton multiplicity (all final-state protons vs protons inside the window).*

## 3. TKI observables

Definitions follow Lu et al., PRC 94 (2016) 015503 and Furmanski & Sobczyk, PRC 95 (2017) 065501, as cited by the paper; ẑ is the neutrino direction, p_T the transverse momenta of the muon (μ) and the leading proton (p):

- δp_T = |p_T^μ + p_T^p|
- δα_T = arccos[ −p̂_T^μ · δp_T / |δp_T| ]
- φ_T = arccos[ −p̂_T^μ · p̂_T^p ]
- δp_Tx = (ẑ × p̂_T^μ) · δp_T, δp_Ty = −p̂_T^μ · δp_T (negative when the proton carries less transverse momentum than the muon)
- δp_L = R/2 − (m_A'² + δp_T²)/(2R), R = m_A + p_L^μ + p_L^p − E^μ − E^p
- p_n = √(δp_T² + δp_L²)

MAT's legacy `MnvRecoShifter::Calc_tki_vars` uses the opposite sign for both δp_Tx and δp_Ty; the paper's Fig. 1 and its released grid match the convention above. Nuclear masses are not stated by the paper; the manifest carries m_A = 11.174864 GeV (¹²C nuclear mass), b = 27.13 MeV (carbon excitation, LE TKI paper arXiv:1805.05486), so m_A' = m_A − m_n + b = 10.262429 GeV, status **default**. The formulas are verified by an exact p_n recovery on a four-momentum-conserving knockout event and by frame-rotation consistency (`tests/test_tki.py`).

TKI of the fiducial truth-level signal (223 595 events; median, 16th–84th percentile):

| observable | median | 16 %–84 % |
|---|---|---|
| δp_T [GeV/c] | 0.264 | 0.091–0.732 |
| δp_Tx [GeV/c] | 0.000 | -0.226–0.227 |
| δp_Ty [GeV/c] | -0.124 | -0.619–0.070 |
| δα_T [deg] | 132.2 | 47.9–167.4 |
| φ_T [deg] | 16.7 | 2.5–87.0 |
| δp_L [GeV/c] | 0.092 | -0.022–0.274 |
| p_n [GeV/c] | 0.349 | 0.126–0.767 |

![truth-level TKI](figs/pl1A_signal_tki.png)
*Figure 2 — TKI distributions of the fiducial truth-level signal: the Fermi peak in p_n near 0.2 GeV/c with the FSI tail, δα_T rising toward 180°, δp_Ty peaked at zero with the deceleration tail, φ_T peaked at zero.*

## 4. Event selection (reconstruction level)

### 4.1 The paper's selection and its transcription

The paper (Sec. "Analysis and results"): events have a negatively charged muon reconstructed in MINOS and at least one proton candidate whose energy is measured by range; protons that exit or interact inelastically are flagged by the end-of-track deposits and stopping protons with a Bragg-peak hit pattern are accepted, the Bragg-peak shape also vetoing pions; no Michel electron candidates near the vertex or any track endpoint; no more than one isolated cluster of energy. The muon and the highest-momentum proton candidate define the TKI variables. A commented-out tex line (unpublished) gives the reconstruction windows: muon 17°, 2–20 GeV/c; proton 90°, 400–1300 MeV/c, "a larger fiducial to avoid cutting on the efficiency edge".

Transcription onto the MasterAnaDev AnaTuple (selection `minerva_ccqelike_1mu1p_v0`, values in `selection.params` of the channel manifest):

| # | cut | paper wording | tuple branch / value | status |
|---|---|---|---|---|
| 1 | ZRange, Apothem | fiducial vertex in the tracker | `vtx` in z 5980–8422 mm, hexagon apothem 850 mm (inclusive channel's box) | open (paper fiducial not stated) |
| 2 | HasMINOSMatch, NoDeadtime, IsNeutrino | μ⁻ reconstructed in MINOS | `isMinosMatchTrack == 1`; ≤ 1 dead discriminator pair upstream; `MasterAnaDev_minos_trk_qp < 0` | decided (inclusive chain) |
| 3 | MuonWindow | 17°, 2–20 GeV/c | beam-frame `muon_thetaX/Y` angle < 17°, 2 < p < 20 GeV/c | default (commented tex line) |
| 4 | HasProtonCandidate | ≥ 1 proton candidate, energy by range | `MasterAnaDev_proton_P_fromdEdx > 0` (the tool's primary candidate; secondary candidates never exceed it) | decided |
| 5 | ProtonContained | exiting protons not well reconstructed | `MasterAnaDev_hadron_isExiting[0] == 0` (slot 0 = proton candidate in ~93 % of events) | default |
| 6 | ProtonScore | Bragg-peak stopping proton, pion veto | `MasterAnaDev_proton_score1 > 0.35` | default (paper gives no value) |
| 7 | ProtonWindow | 90°, 400–1300 MeV/c | `MasterAnaDev_proton_theta` (beam frame) < 90°, 0.4 < P_fromdEdx < 1.3 GeV/c | default (commented tex line) |
| 8 | NoMichel | no Michel electron near vertex or track ends | `improved_nmichel == 0` | decided |
| 9 | IsoBlobs | ≤ 1 isolated cluster | `n_nonvtx_iso_blobs ≤ 1` (`_all` variant agrees in 95 % of candidate events) | decided (counter choice open) |

`MasterAnaDev_proton_theta` equals the angle of the detector-frame (Px, Py, Pz)_fromdEdx after the NuMI beam rotation (median difference 0.0000° on the certification file), the same convention as the muon.

### 4.2 Cutflow, purity, efficiency

Cumulative counts; MC scaled by the POT ratio 0.22028; purity = fraction of selected MC that is truth signal inside the fiducial volume; efficiency = selected truth signal / 223 595 (the fiducial signal count of Sec. 2.2; `cutflow.json` of the selection run):

| step | data | MC (scaled) | data / MC | MC purity | MC efficiency |
|---|---|---|---|---|---|
| all_candidates | 2 791 649 | 1 675 406.0 | 1.666 | 0.024 | 0.802 |
| ZRange | 1 160 098 | 1 135 743.8 | 1.021 | 0.035 | 0.796 |
| Apothem | 712 527 | 719 831.9 | 0.990 | 0.054 | 0.792 |
| HasMINOSMatch | 397 918 | 435 190.6 | 0.914 | 0.087 | 0.772 |
| NoDeadtime | 390 745 | 429 350.4 | 0.910 | 0.088 | 0.763 |
| IsNeutrino | 367 061 | 401 902.2 | 0.913 | 0.092 | 0.748 |
| MuonWindow | 326 877 | 363 045.5 | 0.900 | 0.101 | 0.742 |
| HasProtonCandidate | 183 384 | 204 192.7 | 0.898 | 0.126 | 0.521 |
| ProtonContained | 177 438 | 197 573.4 | 0.898 | 0.129 | 0.518 |
| ProtonScore | 64 924 | 73 705.0 | 0.881 | 0.259 | 0.387 |
| ProtonWindow | 61 411 | 69 489.1 | 0.884 | 0.268 | 0.378 |
| NoMichel | 47 800 | 54 736.5 | 0.873 | 0.328 | 0.365 |
| IsoBlobs | **30 049** | **35 244.6** | **0.853** | **0.484** | **0.347** |

The paper quotes, for the CH tracker, efficiency 28 % and purity 60 %. This selection: efficiency 34.7 %, purity 48.4 %. The excess of data over MC before the fiducial cuts (1.67) is rock-muon and non-fiducial activity that the MC does not simulate; after the fiducial cuts the ratio is 0.99, 0.90 at the proton-candidate requirement, 0.88 after the score cut and 0.853 ± 0.005 (stat) for the selected sample. The MC is unweighted, so these ratios are not normalisation statements.

Composition of the 160 001 selected MC candidates (categories from the reco rows' truth):

| category | n MC | fraction |
|---|---|---|
| signal QE | 34 648 | 0.217 |
| signal 2p2h | 19 327 | 0.121 |
| signal RES (pion absorbed) | 19 957 | 0.125 |
| signal DIS/other | 3 588 | 0.022 |
| background, single π± | 38 316 | 0.239 |
| background, single π⁰ | 18 898 | 0.118 |
| background, multi-pion | 10 389 | 0.065 |
| background, no pion (out of window / fiducial, NC, ν̄, …) | 14 878 | 0.093 |

Proton-score threshold scan with all other cuts fixed:

| `proton_score1` > | data | MC (scaled) | MC purity | MC efficiency |
|---|---|---|---|---|
| none | 59 191 | 68 055.6 | 0.333 | 0.460 |
| 0.2 | 36 248 | 41 501.4 | 0.457 | 0.385 |
| **0.35** | **30 049** | **35 244.6** | **0.484** | **0.347** |
| 0.5 | 24 078 | 29 158.6 | 0.503 | 0.298 |
| 0.6 | 19 749 | 24 739.4 | 0.513 | 0.258 |
| 0.7 | 14 715 | 19 318.5 | 0.524 | 0.206 |
| 0.8 | 8 654 | 12 012.4 | 0.538 | 0.131 |

No threshold reaches the paper's purity: the plateau near 0.54 means the remaining background (dominantly single π± and π⁰ with the pion undetected) is not separated by the dE/dx score. The paper's Bragg-peak / end-of-track criterion, its isolated-cluster counter, or its Michel tagger may differ from the branches used here, or the sideband tuning may absorb part of the gap. This is the main open question of the selection.

### 4.3 Data versus MC

Selected sample, MC stacked by category and scaled to the data POT, data with Poisson errors, ratio panels with the MC-statistics band. Bin edges are the paper's released grids (verified from `anc/tki_release.root`). Events in range and data/MC per grid (`summary.json` of the selection run):

| measurement | data in range | MC scaled in range | data NaN | MC NaN | data/MC first bin | last bin | range over bins |
|---|---|---|---|---|---|---|---|
| muon p [GeV/c] | 30 049 | 35 244.6 | 0 | 0 | 0.75 | 1.06 | 0.75–1.06 |
| muon θ [deg] | 30 049 | 35 244.6 | 0 | 0 | 0.44 | 0.91 | 0.44–0.93 |
| muon p_T [GeV/c] | 30 049 | 35 243.9 | 0 | 0 | 0.37 | 0.81 | 0.37–0.99 |
| leading proton p [GeV/c] | 26 485 | 31 534.7 | 0 | 0 | 0.92 | 0.84 | 0.78–0.92 |
| leading proton θ [deg] | 29 441 | 34 792.6 | 0 | 0 | 0.82 | 1.00 | 0.78–1.00 |
| leading proton p_T [GeV/c] | 30 030 | 35 222.6 | 13 | 71 | 0.87 | 1.19 | 0.82–1.19 |
| δp_T [GeV/c] | 30 020 | 35 211.3 | 13 | 71 | 0.87 | 0.73 | 0.73–0.87 |
| δp_T, fine grid [GeV/c] | 30 020 | 35 211.3 | 13 | 71 | 0.82 | 0.76 | 0.76–0.89 |
| δp_Tx [GeV/c] | 30 036 | 35 229.0 | 13 | 71 | 0.74 | 0.75 | 0.74–0.88 |
| δp_Ty [GeV/c] | 30 036 | 35 229.0 | 13 | 71 | 0.66 | – | 0.56–0.90 |
| δα_T [deg] | 30 036 | 35 229.0 | 13 | 71 | 0.75 | 0.86 | 0.75–0.89 |
| φ_T [deg] | 30 036 | 35 229.0 | 13 | 71 | 0.86 | 0.74 | 0.74–0.92 |
| δp_L [GeV/c] | 30 031 | 35 223.9 | 13 | 71 | 1.01 | – | 0.67–1.08 |
| p_n [GeV/c] | 30 036 | 35 229.0 | 13 | 71 | 0.78 | 0.87 | 0.77–0.98 |

![muon and proton kinematics](figs/pl1A_data_vs_mc_muon_proton_2x2.png)
*Figure 3 — Muon momentum and angle, leading-proton momentum and angle.*

![δp_T (fine grid)](figs/pl1A_data_vs_mc_dpt_fine.png)
*Figure 4 — δp_T (fine grid).*

![δp_T (released grid)](figs/pl1A_data_vs_mc_dpt.png)
*Figure 5 — δp_T (released grid).*

![δα_T](figs/pl1A_data_vs_mc_alpha.png)
*Figure 6 — δα_T.*

![φ_T](figs/pl1A_data_vs_mc_phi.png)
*Figure 7 — φ_T.*

![δp_Tx](figs/pl1A_data_vs_mc_dptx.png)
*Figure 8 — δp_Tx.*

![δp_Ty](figs/pl1A_data_vs_mc_dpty.png)
*Figure 9 — δp_Ty.*

![δp_L](figs/pl1A_data_vs_mc_pl.png)
*Figure 10 — δp_L.*

![p_n](figs/pl1A_data_vs_mc_pn.png)
*Figure 11 — p_n.*

![muon p_T](figs/pl1A_data_vs_mc_muon_pt.png)
*Figure 12 — muon p_T.*

![leading proton p_T](figs/pl1A_data_vs_mc_proton_pt.png)
*Figure 13 — leading proton p_T.*

**Reading of the figures (playlist 1A).** With 470 times the statistics of the single certification
file, the shapes still follow the MC: δp_T and p_n peak where the MC peaks, δα_T rises toward 180°,
φ_T is concentrated at small angles, δp_Tx is symmetric and δp_Ty carries the negative tail. The
data/MC ratio is 0.853 ± 0.005 overall but it is not flat: it rises with the muon angle from 0.44 in
the first degree to 0.91 at 15–17° and with the muon momentum from 0.75 (2–3 GeV/c) to 1.06
(14–20 GeV/c); the muon p_T, which tracks Q² at these angles, goes from 0.37 in its first bin to
0.81 in its last. Against that, δp_T is nearly flat (0.73–0.87 on the released grid, 0.76–0.89 on the
fine grid), δα_T moves from 0.75 to 0.86, and the leading-proton momentum falls from 0.92 at
0.5–0.625 GeV/c to 0.84 above 1 GeV/c. The muon-angle and muon-p_T trends are the signature of the
low-Q² (RPA) suppression and the 2p2h enhancement that the MINERvA tune applies and this unweighted
central-value MC lacks; the transverse-imbalance shapes are less sensitive to those weights, which is
why they stay flat. No normalisation should be read from these ratios before the weight set is decided.

## 5. Caveats

- **Exposure.** 8.969 × 10¹⁹ POT of data (30 049 selected events) against the paper's 10.61 × 10²⁰.
- **MC weights.** The official MC is used as generated: no MINERvA tune (2p2h enhancement, RPA, pion-production retune), no flux constraint, no detector-systematic universes. Ratios are shape statements.
- **Unpublished cut values.** The muon and proton reconstruction windows come from commented-out tex lines; the score threshold and the containment criterion are this analysis's defaults.
- **Fiducial volume and normalisation.** The paper does not state its CH fiducial volume or target count; absolute comparisons need them.
- **Frame.** All angles, truth and reco, are in the NuMI beam frame; the detector-frame numbers differ (the 70° proton edge moves by the 3.4° tilt).
- **Truth skims.** The truth-level numbers come from the per-file skims (CC events in an enlarged tracker box with the final-state list reduced to derived columns, `docs/decisions.md` 2026-09-14); the full truth tables stay on DUNE scratch.

## 6. Open questions (tracked in `docs/open_questions.md`)

1. Purity gap versus the paper (0.48 vs 0.60): which score threshold, containment criterion and isolated-blob counter to adopt.
2. The CH fiducial volume and n_nucleons of the paper.
3. Which released grid is the channel's published measurement (14 exist; the placeholder is the leading-proton momentum).
4. Ratification of the TKI constants (m_A, b) and confirmation of the δp_Tx sign convention against the paper's schematic.
5. Events with a second final-state muon (10⁻⁴ level): signal under the primary-lepton definition, background under MAT's `IsQELike`.
6. The MC weight set before any normalisation is read; the playlist dependence of the MC efficiency (full-FHC report).
7. Which playlists to write into the committed channel manifest (`data.playlists`) and where the full truth archives live.

## 7. Reproduce

```bash
cd /exp/dune/data/users/liangliu/ndp-dev/ndp-platform
# products: grid/README.md (ndp grid plan / submit / status / harvest --no-archive; ndp data merge --beam FHC --playlist <pl> --kind data,mc)
# a copy of channels/minerva_me_ccqelike_1mu1p.yaml with data.playlists = {mc: [1A], data: [1A]} and products_dir set:
python -m ndp signal    --channel <manifest.yaml> --slug signal_minerva_me_ccqelike_1mu1p_FHC1A
python -m ndp selection --channel <manifest.yaml> --slug selection_minerva_me_ccqelike_1mu1p_FHC1A
python report/make_report.py --label "playlist 1A" --tag pl1A --signal runs/2026-09-14_signal_minerva_me_ccqelike_1mu1p_FHC1A --selection runs/2026-09-14_selection_minerva_me_ccqelike_1mu1p_FHC1A --out <this file>
```
