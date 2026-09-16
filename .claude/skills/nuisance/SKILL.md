---
name: nuisance
description: "NUISANCE (generator vs generator vs data comparison framework) and nusystematics (cross-section systematic providers wrapping GENIE Reweight) in this project, built up part by part. Part 1 = installation: what is installed where (GENIE Reweight, nusystematics + systematicstools + fhicl-cpp, NUISANCE), verifying it, how the pixi build works, standalone installs, updating, gotchas. Use when the user says 'install NUISANCE', 'build nuisance', 'is NUISANCE installed', 'nusystematics', 'GENIE Reweight', 'nuisflat', 'PrepareGENIE', 'systematicstools', 'fhicl', or when a task needs NUISANCE or nusystematics. Triggers - nuisance, nuis, nuisflat, nuiscomp, PrepareGENIE, PrepareNuWroEvents, nusystematics, nusyst, systematicstools, GENIE Reweight, GENIE_REWEIGHT, DumpConfiguredTweaksNuSyst, fhicl."
argument-hint: "[install|check|update] [--smoke]"
---

# nuisance — NUISANCE and nusystematics in the NDP platform

Skill root: `<ndp-platform>/.claude/skills/nuisance/` (also reachable as `ndp-dev/.claude/skills/nuisance`,
a symlink). Same layout and rules as the `gibuu` skill. Every fact in `references/` was verified on the
FNAL EAF node and carries the date it was checked; re-verify anything that looks stale before acting on it.

## Parts

| part | covers | file | status |
|---|---|---|---|
| 1 Installation | GENIE Reweight, nusystematics, NUISANCE: what is installed, verify, build, standalone, update, gotchas | `references/install.md` | written 2026-09-13 |
| 2 Preparing inputs | PrepareGENIE / PrepareNuWroEvents / PrepareGiBUU, fluxes, target fractions, `nuis prep` | `references/inputs.md` | not written yet |
| 3 Comparisons and flat trees | `nuisflat`, `nuiscomp`, sample cards, `nuis comp`, plotting | `references/comparisons.md` | not written yet |
| 4 Systematics | nusystematics tool configs, parameter headers, weight dumps, use inside NUISANCE | `references/systematics.md` | not written yet |

A part marked *not written yet* does not exist. Do not improvise it from general NUISANCE knowledge as
if it were verified here: say so, help from the tools' own help text (`nuis <verb> help`, running a
program with no arguments, `$NUISANCE_SRC/README.md`, `$NUSYST_SRC/README.md`) with that caveat, and
offer to write the part.

## Part 1 in short (details and every number: `references/install.md`)

- **Installed (2026-09-13):** GENIE Reweight R-1_04_02 in place at `external/genie/Reweight`;
  nusystematics v02_00_07 (+ systematicstools v02_00_03 + standalone fhicl-cpp 4.18.01) in
  `external/nusystematics/install`; NUISANCE main (86c64b44, 2026-08-27, version string 2.9.9) in
  `external/nuisance/install` with GENIE + GENIEReWeight, NuWro, nusystematics, Prob3plusplus and the
  NuHepMC input handler; NEUT off. All against the in-repo GENIE R-3_06_02, NuWro 25.11.1 and pixi ROOT.
- **Environment:** `activate.sh` (pixi) exports `GENIE_REWEIGHT`, `NUSYST`, `NUSYST_SRC`, `NUISANCE`,
  `NUISANCE_SRC`, `nusystematics_ROOT`, `FHICL_FILE_PATH`, `GENIE_XSEC_TUNE`, and puts every `bin` on PATH.
- **Verify:** `bash .claude/skills/nuisance/scripts/nuisance_check.sh [--smoke]` — installs, programs,
  shared-library resolution (including which NuHepMC cpputils copy is loaded), the env conventions, and
  with `--smoke` the platform's `pixi run test-nuisance` (PrepareGENIE / PrepareNuWroEvents → nuisflat,
  nusystematics ResIso weight dump).
- **Build:** `pixi run build-genie-reweight`, `pixi run build-nusystematics`, `pixi run build-nuisance`
  (or `pixi run build-nuisance-stack`); scripts in `scripts/build_*.sh`, idempotent by existence.
- **Origin of the recipe:** the user's `github.com/LiangLiu212/BuildEventGenerators` (FNAL gpvm / UPS).
  UPS `setup` does not work on the EAF node, hence the pixi route; section 1.3 lists every deviation.
- **Gotchas (all handled in the scripts):** gcc 15 vs fhicl-cpp (`-include cassert`), Boost/Eigen/TBB
  added to pixi, ACHILLES's older NuHepMC cpputils shadowing NUISANCE's (same soname), `GENIE_XSEC_TUNE`
  and `.` on `FHICL_FILE_PATH` required by the nusystematics tools, Prepare* tools refuse to overwrite.

## Conventions for growing this skill

- One `references/<part>.md` per part; this file stays a router (parts table + short summaries). When a
  part is added: fill its row above, extend the `description` triggers, keep numbering.
- A reference states only measured or verified values, each with its date and source (command, file,
  URL). Unknowns that the user must supply are written `xxx`, never a plausible-looking value (global
  CLAUDE.md rule).
- Project mechanics live in the repo (`scripts/build_genie_reweight.sh`, `build_nusystematics.sh`,
  `build_nuisance.sh`, `test_nuisance.sh`, `pixi.toml`, `activate.sh`, README section "NUISANCE and
  nusystematics"); references point at them instead of copying them, and get updated when they change.
- `scripts/` here are bash or stdlib-Python only, default to the platform install, and re-exec
  themselves under `pixi run` when the platform environment is not active.
- Gotchas found while working go into the part they belong to, in its final "Gotchas" section.
