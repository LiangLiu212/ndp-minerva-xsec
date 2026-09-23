# 1mu + leading proton: results

Channel `minerva_me_ccqelike_1mu1p` (MINERvA ME FHC, arXiv:2503.15047 signal definition). Model: GiBUU 2025
numu CC on 12C under the NuMI ME FHC flux (energy-scan sample `runs/_generator_cache/gibuu_1f8bcb950be1dcd0`).
Each section adds one stage; the scripts next to this file make the figures.

## 1. GiBUU truth: muon momentum and angle, signal definition only

The only requirement is the truth signal definition of the channel: one muon with angle to the beam below 17 degrees
and momentum between 2 and 20 GeV/c, at least one proton with angle below 70 degrees and momentum between 0.5 and
1.1 GeV/c, no mesons, no baryons heavier than the neutron, no photons above 10 MeV. No reconstruction cut,
efficiency or smearing is applied. Angles are with respect to the beam. Weights are the GiBUU event weights scaled
to the data exposure (1.057e21 POT), so the histograms count expected signal events in the fiducial volume.

Of 979,395 generated CC events, 56,132 pass the signal definition, corresponding to 533,040 signal events at the
data exposure out of 8,902,466 CC events.

![muon momentum, truth signal](figs/muon_p_truth_signal.png)

True muon momentum of the signal events: median 5.25 GeV/c, mean 5.43 GeV/c, peak bin 5.0 to 5.5 GeV/c.

![muon angle, truth signal](figs/muon_theta_truth_signal.png)

True muon angle to the beam of the signal events: median 7.25 degrees, mean 7.67 degrees, peak bin 6.0 to 6.5 degrees.

![muon cos theta, truth signal](figs/muon_costheta_truth_signal.png)

The same angle as cos(theta), from the signal window's edge cos(17 degrees) = 0.9563 to 1 in 44 bins of 0.001: median 0.9926, mean 0.9892, peak bin 0.9970 to 0.9980.

Script: `make_muon_truth_signal.py` (default pixi environment); histogram contents and the numbers above:
`muon_truth_signal.json`.

## 2. Selection cuts applied

The cuts are the channel's reconstruction-level selection (`minerva_ccqelike_1mu1p_v0`): tracker fiducial vertex,
MINOS-matched negative muon, dead time, muon window (17 degrees, 2 to 20 GeV/c), a contained proton candidate with
dE/dx score above 0.35 inside 90 degrees and 0.4 to 1.3 GeV/c, no Michel electron, at most one isolated blob. A
generator sample has no reconstruction, so the cuts are applied as the selection efficiency learned from the official
MC (the 12 ME FHC playlists): per true cell of a grid, efficiency = selected truth signal / truth signal in the
fiducial volume. Each signal event of section 1 is weighted by the efficiency of its true cell on the grid of the
plotted variable; the lower panel of each plot shows the resulting efficiency per bin. Events are otherwise unchanged,
so the distributions are still in true kinematics.

The cuts keep 175,118 of the 533,040 signal events at the data exposure, a mean efficiency of 0.33.

![muon momentum, cuts applied](figs/muon_p_cuts_applied.png)

Muon momentum: the efficiency rises from 0.24 in the 2 to 3 GeV/c bin to 0.38 at 7.5 to 10 GeV/c and falls back to
0.35 above 14 GeV/c; after the cuts the median is 5.25 GeV/c and the mean 5.65 GeV/c, peak bin 5.0 to 5.5 GeV/c.

![muon angle, cuts applied](figs/muon_theta_cuts_applied.png)

Muon angle: the efficiency falls steadily with angle, from 0.46 below 1 degree to 0.13 in the 14 to 17 degree bin,
because the MINOS match favours forward muons; after the cuts the median is 6.25 degrees and the mean 6.79 degrees,
peak bin 5.5 to 6.0 degrees.

![muon cos theta, cuts applied](figs/muon_costheta_cuts_applied.png)

The same in cos(theta): efficiency from 0.09 at the window edge to 0.40 in the most forward bin; after the cuts the
median is 0.9935 and the mean 0.9915, peak bin 0.9970 to 0.9980.

Script: `make_muon_cuts_applied.py` (default pixi environment); histogram contents, the grid efficiencies and the numbers
above: `muon_cuts_applied.json`.

## 3. VBLL surrogate smearing applied

The detector resolution is applied with the ported VBLL surrogate `x60_het` (a VBLL_SurrogateModel checkpoint trained
on MINERvA StandardMC pairs of true and reconstructed muon and proton 4-vectors; port and certification in
`report/VBLL_surrogate_1mu1p.md`). For every signal event of section 2 the true muon and leading-proton 4-vectors
are rotated into the model's detector frame and 20 reconstructed copies are drawn from the model's predictive
distribution; the copies are conditioned on the selection's kinematic windows (85.6 percent of them pass, the
event's weight being shared among the passing ones, no event lost) and each event keeps the selection efficiency of
section 2. The plotted quantities are now the smeared, reconstructed muon momentum and angle to the beam. This is
the prediction that the platform compares with the data in reconstructed space; the lower panel shows the ratio of
the smeared to the unsmeared distribution, the effect of the migration alone.

The smearing conserves the totals of section 2 (175,118 events at the data exposure on the momentum grid).

![muon momentum, VBLL smeared](figs/muon_p_vbll_smeared.png)

Muon momentum: median 5.25 GeV/c and peak bin 5.0 to 5.5 GeV/c as before the smearing; the migration is within 5
percent below 6 GeV/c, removes 5 to 15 percent between 6 and 8 GeV/c, adds up to 70 percent between 9 and 11 GeV/c and
removes 20 to 60 percent above 16 GeV/c. This is the effect of the model's muon energy width, 1.35 GeV on average,
which is larger than the 0.48 GeV core width of the residual it was trained on.

![muon angle, VBLL smeared](figs/muon_theta_vbll_smeared.png)

Muon angle: median 6.25 degrees, peak bin 5.0 to 5.5 degrees (5.5 to 6.0 before); both ends of the window gain,
30 to 70 percent below 2 degrees and 60 to 100 percent above 14 degrees, while 4 to 12 degrees lose about 10 percent.

![muon cos theta, VBLL smeared](figs/muon_costheta_vbll_smeared.png)

The same in cos(theta): median 0.9945 (0.9935 before), peak bin 0.9970 to 0.9980 unchanged; the bins nearest the
window edge, cos(theta) below 0.962, gain a factor 2 to 3 and the most forward bin gains 37 percent.

### Migration matrices

The same smeared copies, histogrammed in two dimensions: true value on the horizontal axis, VBLL-reconstructed value
on the vertical axis. Each column is normalised to the events of its true bin, so a cell reads as the probability
P(reconstructed bin | true bin) in percent. The left panel uses the analysis grid of the corresponding measurement,
the grid on which the platform folds; the right panel uses the fine bins of the plots above. No copy leaves the
grids' ranges, because for these three variables the ranges coincide with the reconstruction windows the copies are
conditioned on.

![muon momentum, VBLL migration](figs/muon_p_vbll_migration.png)

Muon momentum on the 8-bin grid: 32 to 35 percent of the events stay in their true bin below 6 GeV/c, 41 percent at
6 to 7.5 GeV/c and 51 to 57 percent in the three wide bins above 7.5 GeV/c; averaged over the signal, 37 percent stay.

![muon angle, VBLL migration](figs/muon_theta_vbll_migration.png)

Muon angle on the 14-bin grid (1-degree bins up to 12 degrees): 40 to 45 percent stay in their bin below 3 degrees,
falling to 10 percent at 10 to 12 degrees where the smearing exceeds the bin width, 17 to 20 percent in the two wider
bins above 12 degrees; 23 percent on average.

![muon cos theta, VBLL migration](figs/muon_costheta_vbll_migration.png)

Muon cos(theta) on the 8-bin grid: 80 percent stay in the most forward bin (cos(theta) above 0.9945) and 7 to 23
percent in the six bins below 0.989; 47 percent on average.

The matrices in percent, the diagonal fractions and the per-column losses are in `muon_vbll_smeared.json` under
`migration`.

### Migration matrices from the official MC, same bins

The same matrices from the MINERvA StandardMC of the 12 ME FHC playlists, the sample the platform's binned responses
are learned from: reconstructed candidates passing the channel selection whose truth is signal inside the fiducial
volume, 891,177 events, true muon kinematics from the primary lepton and reconstructed ones from the MINERvA
reconstruction (`reco_p`, `reco_theta`). Same analysis grids, same fine bins, same column normalisation. On the analysis
grids the matrices reproduce the migration counts stored in the binned responses exactly.

![muon momentum, MC migration](figs/muon_p_mc_migration.png)

![muon angle, MC migration](figs/muon_theta_mc_migration.png)

![muon cos theta, MC migration](figs/muon_costheta_mc_migration.png)

![diagonal fractions, MC vs VBLL](figs/migration_diagonal_mc_vs_vbll.png)

Probability of reconstructing in the true bin, official MC against the VBLL surrogate, averaged over the signal:

| grid | official MC | VBLL x60_het |
|---|---|---|
| muon momentum, 8 bins | 68 % (63 to 83 % per bin) | 37 % (32 to 57 %) |
| muon angle, 14 bins | 77 % (72 to 92 %) | 23 % (10 to 45 %) |
| muon cos(theta), 8 bins | 91 % (71 to 96 %) | 47 % (7 to 80 %) |

The MINERvA reconstruction keeps three quarters of the muons inside a 1-degree bin and two thirds inside a 1 GeV/c
bin; the VBLL surrogate, whose sigma is calibrated to the root-mean-square of a heavy-tailed residual rather than to
its core, spreads them over several bins, most strongly in angle. Matrices, counts and diagonals: `muon_mc_migration.json`;
script `make_muon_mc_migration.py` (default pixi environment, one playlist at a time, about 4 minutes).

Script: `make_muon_vbll_smeared.py` (ml pixi environment: `pixi run -e ml python ...`); histogram contents, including
the smeared distributions without the efficiency, and the numbers above: `muon_vbll_smeared.json`.
