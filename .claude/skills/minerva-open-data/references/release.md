# Part 1 — The MINERvA Open Data Product and how to get it

Verified 2026-09-13 from https://minerva.fnal.gov/opendata/ and its sub-pages, `xrdfs` listings of the
XRootD area, the published file lists, and the MinervaExpt GitHub repositories. The page announces
further tuple fixes "later in 2026" (flux information in the tuples, EM-shower momentum bias,
ML vertex consistency) and Low Energy data "soon": re-check before quoting numbers.

## 1.1 What it is

- "All the neutrino and antineutrino data from the MINERvA experiment" in two beam eras: Low Energy
  (LE, flux peak ≈ 3–3.5 GeV, 2009–2012) and Medium Energy (ME, peak ≈ 6 GeV, 2013–2019), each in
  Forward Horn Current (FHC, neutrino-dominated) and Reverse Horn Current (RHC, antineutrino) modes.
  As of 2026-09-13 the ME FHC and ME RHC playlists are released; the `LowEnergy_FHC/{Data,MC}` and
  `LowEnergy_RHC` directories exist on XRootD but are empty.
- "The main data release are .root files for data and simulation", "relatively flat tuples … that
  contain a variety of reconstructed objects", produced by MINERvA's `MasterAnaDev` (MAD) tool. The
  data "has already been pre-selected to contain either a muon candidate or an electron candidate".
- Simulation comes with "a way to access our uncertainties on that simulation, including flux,
  neutrino interaction, and detector uncertainties" (GENIE weight branches in the tuples, flux and
  reweight files, and the MAT software).
- DOI `10.15484/3022562` (resolves to the opendata page). License CC0. Citation: the DOI (BibTeX at
  https://inspirehep.net/literature/3128404) plus the detector NIM paper (Nucl. Inst. Meth. A743 (2014)
  130, arXiv:1305.5199); suggested acknowledgment: "The authors thank the MINERvA Collaboration for
  making their data, simulated data, and analysis tools available to the community." "We do not
  require MINERvA review of non-MINERvA manuscripts"; support is "limited" and not guaranteed;
  contact the spokespeople (https://minerva.fnal.gov/collaboration/).

## 1.2 Playlists (run periods)

A playlist is a data-taking period with fixed beam and detector conditions (https://minerva.fnal.gov/minerva-playlists/).
Rules stated there: process all playlists of ME FHC or all of ME RHC, "never mix the two"; LE and
ME "can't be mixed because the flux and other conditions are too different"; each playlist gets "a
flux weight suitable for the entire run period". The helium-target column is the *simulated* status;
in data the target was transitioning during 1B, 1O, 6G.

| ME FHC playlist | water target | helium target (MC) | data POT (×10²⁰) |
|---|---|---|---|
| minervame1A | empty | empty | 0.90 |
| minervame1B | empty | full* | 0.19 |
| minervame1C | empty | full | 0.43 |
| minervame1D | empty | full | 1.4 |
| minervame1E | empty | full | 1.0 |
| minervame1F | empty | full | 1.7 |
| minervame1G | empty | empty | 1.4 |
| minervame1L | full | empty | 0.13 |
| minervame1M | full | empty | 2.1 |
| minervame1N | full | empty | 1.1 |
| minervame1O | full | empty* | 0.30 |
| minervame1P | full | full | 0.47 |

| ME RHC playlist | water target | helium target (MC) | data POT (×10²⁰) |
|---|---|---|---|
| minervame5A | full | empty | 0.55 |
| minervame6A | full | full | 1.6 |
| minervame6B | full | full | 1.0 |
| minervame6C | full | full | 1.1 |
| minervame6D | full | full | 1.2 |
| minervame6E | full | full | 0.87 |
| minervame6F** | full | full | 1.4 |
| minervame6G | full | empty* | 0.74 |
| minervame6H | full | empty | 1.1 |
| minervame6I | full then empty | empty | 0.79 |
| minervame6J | empty | empty | 0.91 |

`*` transitioning in data; `**` "A bug in the Playlist 6F Standard Monte Carlo sample machine learning
vertex predictions was found Apr. 2, 2026" — fix pending. Playlist 6I is additionally split into
`6I_WaterFull` / `6I_WaterEmpty` file lists (same files, categorised by water status).

## 1.3 Files: naming, layout, sizes

XRootD endpoint (anonymous; three sanctioned ways: stream with xrootd, copy with `xrdcp`, or mount
PNFS inside Fermilab — "Direct use is not advised as it will overload the cache (dCache) system"):

```
root://fndcadoor.fnal.gov:1095//pnfs/fnal.gov/usr/minerva/persistent/OpenData/<path>
```

Layout seen with `xrdfs fndcadoor.fnal.gov:1095 ls` on 2026-09-13:

```
OpenData/FluxAndReweightFiles/      README.txt + FluxAndReweightFiles_Tarred_Feb_20_2026_1145_FNALTime.tgz (2.6 GB)
OpenData/LowEnergy_FHC/{Data,MC}    empty (announced "soon");  OpenData/LowEnergy_RHC/
OpenData/MediumEnergy_FHC/Data/Playlist1A … Playlist1P/
OpenData/MediumEnergy_FHC/MC/{StandardMC,Extended2p2h,CCDiffractivePion,NCDiffractivePion,CCCoherentPion,
                              DSCalorimeters,ElectronNeutrino,NeutrinoElectronElastic}/Playlist<X>/
OpenData/MediumEnergy_RHC/Data/…, MediumEnergy_RHC/MC/<same eight sample dirs>/
```

File names: `MasterAnaDev_data_AnaTuple_run<NNNNNNNN>_Playlist.root` and
`MasterAnaDev_mc_AnaTuple_run<NNNNNNNN>_Playlist.root`. Data run numbers are the real ones
(Playlist 1A: 253 files, runs 6038–10066, 26 MB–570 MB each); MC "runs" are generation jobs
(FHC 1A StandardMC: 41 files, runs 110000–110040, ~21.5 GB each; RHC 6A StandardMC starts at run
122000, 64 files). Special samples reuse the StandardMC run numbering inside their own directory
(e.g. `MC/Extended2p2h/Playlist1A/…run00110000…`, 4 files for 1A; `MC/ElectronNeutrino/Playlist1A`, 41).
Everything published is 4584 list entries (see below); the MINERvA-101 wiki quotes "approximately 106 GB
total (2 GB data, 102 GB MC files, 4 GB supporting materials)" for its tutorial subset.

File lists: one plain-text file per (beam, sample, playlist) with one xrootd URL per line, hosted on
the site (`https://minerva.fnal.gov/wp-content/uploads/2026/03/MediumEnergy_<FHC|RHC>_<Data|StandardMC>_Playlist<X>.txt`;
special samples under `uploads/2025/10/MediumEnergy_FHC_<Extended2p2h|CCDiffractivePion|DSCalorimeters|
NeutrinoElectronElastic|ElectronNeutrino>_Playlist<X>.txt`, note `DSCalorimeters_Playlist1A-1.txt`).
116 lists on 2026-09-13 with 4584 entries (the two 6I water-split lists repeat 6I's files); the
standard lists name the door `root://fndcadoor.fnal.gov:1095/…`, the special-sample lists
`root://fndca1.fnal.gov:1095/…` (both serve the same namespace). Per kind: FHC Data 1818, FHC
StandardMC 489, RHC Data 1085, RHC StandardMC 473, ElectronNeutrino 400, DSCalorimeters 236,
Extended2p2h 44 (≈25.3 GB each), CC/NC DiffractivePion 13 each, NeutrinoElectronElastic 13.
`scripts/minerva_od.py filelists` scrapes and downloads them; the MINERvA-101
`runEventLoop <data.txt> <mc.txt>` consumes exactly this format.

## 1.4 Standard and special Monte Carlo

- **StandardMC**: the official simulation matching each playlist (POT, target status). Generator tag
  is *not* stored in the tuples (`Meta` has only POT and entry counts); GENIE 2.12.6 is asserted in
  this workspace's channel YAML but marked unconfirmed in both repos' open questions.
- **Extended 2p2h** (`MC/Extended2p2h`): Valencia 2p2h "modified … to produce cross section and events
  up to q3 = 2 GeV", "roughly 16 times the predicted data statistics". Use: "configure your code to
  weight the other 2p2h sample to zero, and to weight this sample by the ratio of generated POT";
  weight the unphysical Q²≈0 high-ω region "down to zero" when applying original-Valencia predictions.
- **Coherent and diffractive pion** (`MC/CCDiffractivePion`, `MC/NCDiffractivePion`; `MC/CCCoherentPion`
  exists on XRootD, the page says CC coherent "will be available in the future"): GENIE's Rein
  diffractive model (off by default in GENIE), "roughly 30 times the predicted data statistics".
- **Downstream calorimeters** (`MC/DSCalorimeters`): events originating in ECAL/HCAL with the tracker
  reconstruction, ~4× data; "use this sample OR the regular sample" (overlap double/undercounts planes
  at the end of the tracker); HCAL analyses untested, calibration care needed.
- **Neutrino–electron elastic** (`MC/NeutrinoElectronElastic`): the flux standard candle, ~300× data.
- **Electron neutrino** (`MC/ElectronNeutrino`): νe are ~1 % of the flux; ~4× data; not every playlist
  has a list (1F missing on the page).

## 1.5 Flux, reweight and parameter files

- https://minerva.fnal.gov/minerva-fluxes/: ME (12E20 POT each mode) recommended flux from the
  simultaneous ν–e + inverse-muon-decay fit (Zazueta et al., PRD 107 (2023)):
  `https://minerva.fnal.gov/wp-content/uploads/2023/06/improved_flux_constrain_datarelease_excelformat.zip`
  (νμ, ν̄μ, νe, ν̄e fluxes, uncertainties, covariances, README); nuclear-target fluxes
  `…/2023/06/MINERvA_ME_NuMU_Flux_Nuclear_Targets_Constrained.root_.zip` (per-target, ±10 % on the
  falling edge, PRL 130, 161801); LE fluxes (PRD 94, 092005 unconstrained; ν-mode constrained only).
  The 2019 ν-mode-only constrained flux (PRD 100, 092001) is the "legacy" one. This workspace uses the
  arXiv:2110.13372 supplemental Table I flux (νμ/m²/10⁵ POT/GeV, 0–100 GeV) — see part 3.
- `OpenData/FluxAndReweightFiles/FluxAndReweightFiles_Tarred_Feb_20_2026_1145_FNALTime.tgz` (2.6 GB;
  README: "the parameter and reweight files used by MINERvA Analyses", dated tarballs replace earlier
  ones). Its first entries are `MATFluxAndReweightFiles/{flux/flux-g4numiv5-pdg{14,-14}-minerva<1,5,13,LE-FHC>.root,
  flux/sys/*_MnvH1D.root, dataMCWiggle/…}`; the rest of the archive (MParamFiles are expected there by
  the MINERvA-101 wiki) was not listed — inspect with `tar tzf` after downloading. These are MAT
  `MnvH1D` files: unreadable without the MAT libraries (part 3).
- Tunes (https://minerva.fnal.gov/minerva-tunes/): MnvTune vX.Y.Z, "not cumulative"; v1 = Valencia RPA
  + low-recoil 2p2h fit, v2 adds the low-Q² pion suppression, v3 SuSA 2p2h + Bodek-Ritchie tail,
  v4 full bubble-chamber pion fit (CCNormRes 1.15, MaRES 0.94); Y variants for coherent norm, pion
  suppression, axial form factor, nuclear model, KE weights. Applied through MAT reweighters, not
  stored in the tuples.

## 1.6 Official software (github.com/MinervaExpt)

| repo | role |
|---|---|
| `MAT` | MINERvA Analysis Toolkit (arXiv:2103.08677): many-universe systematics, `MnvH1D/MnvH2D` |
| `MAT-MINERvA` | MINERvA plugins: `FluxReweighter`, `MuonFunctions`, `CCInclusiveCuts`, standard systematics; installs MAT + UnfoldUtils |
| `UnfoldUtils` | RooUnfold fork with systematics propagation |
| `GENIEXSecExtract` | closure tests for event loops |
| `MINERvA-101-Cross-Section` | the tutorial: `runEventLoop <data.txt> <mc.txt>` → `ExtractCrossSection`; basis of PRD 104, 092007 (ME FHC inclusive) |
| `Tuple-Documentation` | `MAD_tuple_MainDoc.csv/.xlsx`, 257 documented branches (groups: part, MasterAnaDev, truth, gamma1/2, pi0, shower, sec, event, hadron, proton, slice…) — "work in progress" |
| `CCQENu`, `CC-CH-pip-ana`, `NucCCNeutrons`, `LowRecoilPions`, `NuETKI`, … | analysis packages on the same tuples |

MINERvA-101 wiki prerequisites: ROOT 6 (6.28.12 tested, 6.38 works) built with XRootD, `xrootd-client`
for `xrdcp`, CMake, CVS (for MParamFiles), then MAT, MAT-MINERvA, UnfoldUtils, GENIEXSecExtract,
MParamFiles, MATFluxAndReweightFiles (`xrdcp root://fndca1.fnal.gov:1095/pnfs/fnal.gov/usr/minerva/persistent/OpenData/FluxAndReweightFiles/FluxAndReweightFiles_Tarred_Feb_20_2026_1145_FNALTime.tgz ./ && tar -zxvf …`).
None of this C++ stack is built in this workspace (part 3 explains the Python route used instead).
Also linked: Arachne event display (Fermilab-internal), GENIE2 generation page
(https://minerva.fnal.gov/genie2/), NuMI beam paper (arXiv:1507.06690), test-beam calibration
(arXiv:1501.06431).
