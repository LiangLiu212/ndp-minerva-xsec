---
name: minerva-open-data
description: "MINERvA Open Data Product (minerva.fnal.gov/opendata): what the release contains (Medium Energy FHC/RHC playlists, data + StandardMC + special MC samples, flux and reweight files), how to get it over XRootD, the MasterAnaDev AnaTuple format (trees, POT, key branches, truth, GENIE weight universes, frames, recoil energies), and how data and MC are processed in this workspace (ndp-platform caches/channels/surrogates and the exploration repo's tools) versus the official MAT tutorial. Use when the user says 'MINERvA open data', 'AnaTuple', 'MasterAnaDev', 'playlist', 'minervame1A', 'download MINERvA data', 'POT', 'MAT', 'MnvH1D', 'flux files', 'process the MINERvA MC', or when a task reads or caches MINERvA files. Triggers - minerva, MINERvA, open data, opendata, AnaTuple, MasterAnaDev, MAD, playlist, minervame, fndcadoor, xrdcp, POT_Used, Truth tree, MAT, PlotUtils, MnvH1D, FluxAndReweightFiles, MParamFiles."
argument-hint: "[release|format|processing] [--playlist 1A]"
---

# minerva-open-data — MINERvA's data-preservation AnaTuples, and how we process them

Skill root: `<ndp-platform>/.claude/skills/minerva-open-data/` (also `ndp-dev/.claude/skills/minerva-open-data`,
a symlink). Same layout and rules as the `gibuu` and `nuisance` skills. Every fact in `references/`
was verified on 2026-09-13 (web pages, XRootD listings, the local files) unless dated otherwise;
the release is being revised during 2026, so re-check the page before quoting it.

## Parts

| part | covers | file | status |
|---|---|---|---|
| 1 The release and access | what the Open Data Product is, playlists and POT, standard and special MC, file naming and XRootD layout, file lists, flux and reweight files, tunes, official software, citation | `references/release.md` | written 2026-09-13 |
| 2 Tuple format | trees, Meta/POT, key reco and truth branches, selection recipe, signal definition, GENIE weight universes, units, frames, recoil energies, Truth-cycle gotcha, data/MC scaling | `references/format.md` | written 2026-09-13 |
| 3 Processing here | the two routes (official MAT C++ vs this workspace's Python), download, ndp-platform commands and caches, certified numbers, normalisation, what is not reproducible, adding a playlist | `references/processing.md` | written 2026-09-13 |
| 4 Systematics with MAT | building MAT/MAT-MINERvA/UnfoldUtils, MnvH1D universes, MParamFiles, flux universes, running MINERvA-101 end to end | `references/mat.md` | not written yet |

A part marked *not written yet* does not exist. Do not improvise it from general MINERvA knowledge as
if it were verified here; say so, point at the official sources listed in part 1, and offer to write it.

## Quick answers

- **Where the data is:** `root://fndcadoor.fnal.gov:1095//pnfs/fnal.gov/usr/minerva/persistent/OpenData/`
  (anonymous XRootD, ~85 MB/s from the EAF node). Layout `MediumEnergy_{FHC,RHC}/{Data,MC/<sample>}/Playlist<X>/
  MasterAnaDev_{data,mc}_AnaTuple_run<8 digits>_Playlist.root`; the site publishes one file list per
  playlist and sample (116 lists, 4584 entries on 2026-09-13). One StandardMC file is ~21.5 GB.
- **What is on disk here:** `/exp/dune/data/users/liangliu/minerva-open-data/` holds one ME FHC 1A pair
  (data run 10066, 188 MB; StandardMC run 110040, 21.6 GB) plus `cache/*.npz` built by ndp-platform.
- **Helper:** `python3 .claude/skills/minerva-open-data/scripts/minerva_od.py {filelists,fetch,inspect,status}`
  — list/download the published file lists, xrdcp files with size checks, print a tuple's trees/POT,
  show what is local.
- **Processing in this workspace:** `cd <ndp-platform> && python3 -m ndp data status` /
  `python3 -m ndp data cache --channel minerva_me_cc_inclusive_ptpz` (uproot adapter,
  `ndp/adapters/minerva_anatuple.py`); the exploration repo `ndp-minerva-data-release-exploration`
  holds the PyROOT tools, the data dictionary and the audited cross-section runs. The official route
  (MAT + MINERvA-101 tutorial, C++) is not built here.
- **Three gotchas to keep in mind:** reco muon angles are beam-frame while `mc_primFSLepton` is
  detector-frame (rotate by −0.05887 rad about x); the MC `Truth` tree has two key cycles (use the
  highest); the seven recoil-energy branches fall in three families and none is yet declared canonical.

## Conventions for growing this skill

- One `references/<part>.md` per part; this file stays a router. Add a part: fill its row, extend the
  `description` triggers, keep numbering.
- Only measured or verified values, each dated and sourced (URL, command, file:line). Unknowns the user
  must supply are `xxx`, never a plausible-looking value (global CLAUDE.md rule).
- Project mechanics stay in the repos (`ndp-platform/ndp/adapters/minerva_anatuple.py`,
  `channels/*.yaml`, `tests/test_minerva_certification.py`, the exploration repo's `docs/` and
  `scripts/download_data.sh`); references point at them and get updated when they change. Never edit
  the exploration repo's audited runs from here.
- `scripts/` are stdlib Python (uproot optional, present in the default and platform pixi pythons).
- Gotchas found while working go into the part they belong to, in its final "Gotchas" section.
