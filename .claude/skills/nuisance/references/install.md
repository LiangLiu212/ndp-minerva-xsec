# Part 1 — Installing NUISANCE and nusystematics

Verified on 2026-09-13 on the FNAL EAF node (`jupyter-liangliu-*`, EL9, see the `fnal-eaf-node` memory)
unless another date is given. Sources: the build scripts and logs under `ndp-platform/scripts/` and
`external/`, the projects' CMake files, `nuis config`, and the user's recipe repository
`https://github.com/LiangLiu212/BuildEventGenerators` (fork of afropapp13's; FNAL gpvm + UPS).

## 1.0 What the two packages are

- **NUISANCE** (`github.com/NUISANCEMC/nuisance`, JINST 12 P01016): compares neutrino event generators
  with each other and with published cross-section data. Generator output is first converted into
  NUISANCE's common event format (`PrepareGENIE`, `PrepareNuWroEvents`, `PrepareGiBUU`, `nuis prep`),
  then used for generator-agnostic flat trees (`nuisflat`), data/MC comparisons (`nuiscomp`, `nuis comp`),
  fits and systematic studies (`nuismin`, `nuissyst`) over the built-in experimental samples (ANL,
  ArgoNeuT, BEBC, BNL, electron scattering, FNAL, GGM, K2K, MicroBooNE, MINERvA, MiniBooNE, SciBooNE, T2K).
- **nusystematics** (`github.com/NuSystematics/nusystematics`): cross-section systematics for GENIE 3
  events as "systematic providers" in the `systematicstools` framework (FHiCL-configured; the DUNE/SBN
  convention). The `GENIEReWeight` provider wraps GENIE Reweight; others (shipped in v02_00_07) are
  BeRPAWeight, CCQERPAReweight, DIRT2_Emiss, EbLepMomShift, FSILikeEAvailSmearing, FSIReweight,
  MINERvAE2p2h, MINERvAq0q3Weighting, MiscInteractionSysts, MKSinglePiTemplate,
  NOvAStyleNonResPionNorm, ResIso, ZExpPCAWeighter. Two-step configuration: a *tool configuration*
  FHiCL → `GenerateSystProviderConfigNuSyst` → *parameter headers* FHiCL → `DumpConfiguredTweaksNuSyst`
  (or a client such as NUISANCE) computes per-event response weights.
- **GENIE Reweight** (`github.com/GENIE-MC/Reweight`) is the prerequisite of both: the GENIE weight
  calculators (`libGRwFwk`, `libGRwIO`, `libGRwClc`) live outside the Generator since GENIE 3.

## 1.1 What this project already has

Everything is under `ndp-platform/external/` (gitignored, rebuildable from `scripts/`), built in the
pixi environment against the in-repo GENIE Generator R-3_06_02, NuWro 25.11.1 and pixi ROOT 6.40.04
with gcc 15.2.0 and CMake 4.4.3.

| component | version | where | size | build |
|---|---|---|---|---|
| GENIE Reweight | 1.04.02 (tag R-1_04_02, 2025-06-30) | `external/genie/Reweight` (in place: `lib/`, `bin/`, `src/`) | 4.1 MB | `pixi run build-genie-reweight`, 44 s with `-j16` |
| nusystematics | tag v02_00_07 (2026-01-12; its CMake still says `VERSION 02.00.06`, so `setup.nusystematics.sh` and `nuisance_check.sh` report 02.00.06 — upstream did not bump it); pins systematicstools v02_00_03 and the fhicl-cpp suite FHICLCPP_SUITE_v4_18_01 (cetmodules 3.22.01, cetlib-except, hep-concurrency, cetlib, fhicl-cpp) | source `external/nusystematics/nusystematics`, install `external/nusystematics/install` (28 programs, 13 libraries; the fhicl tools share the prefix) | 643 MB | `pixi run build-nusystematics`, 169 s for the pass that succeeded |
| NUISANCE | main @ 86c64b44 (2026-08-27), `nuis config --version` = 2.9.9, `git describe` = PR/modern_cmake-343-g86c64b44 | source `external/nuisance/nuisance`, install `external/nuisance/install` (31 programs) | 1.1 GB | `pixi run build-nuisance`, 418 s with `-j16` |

NUISANCE features (`nuis config --features`): nusystematics, GENIE, GENIEReWeight, GENIE3, NuWro,
Prob3plusplus; generators (`--generators`): GENIE, NuWro. The NuHepMC input handler is compiled in as
well (`libInputHandler.so` links `libnuhepmc_cpputils.so`) although `--features` does not list it;
GiBUU input support is always compiled (it reads GiBUU's RootTuple ROOT output). NEUT is off (not
installed). GENIE3 reweighting incl. XSecMEC is on (`OPTIONS: GENIEReWeight: TRUE, XSecMECReWeight: TRUE`
in `external/nuisance/nuisance/configure.log`).

Environment (`activate.sh`, sourced by `pixi run` / `pixi shell`):

| variable | value |
|---|---|
| `GENIE_REWEIGHT`, `GENIE_REWEIGHT_VERSION` | `$NDP_EXTERNAL/genie/Reweight`, `R-1_04_02` |
| `NUSYST`, `NUSYST_SRC`, `NUSYST_VERSION` | `$NDP_EXTERNAL/nusystematics/install`, `.../nusystematics`, `v02_00_07` |
| `nusystematics_ROOT`, `systematicstools_ROOT`, `fhicl_cpp_standalone_ROOT` | all `$NUSYST` (CMake `find_package` prefixes for NUISANCE and other clients) |
| `FHICL_FILE_PATH` | `.:$NUSYST/fcl/fcl:$NUSYST_SRC/fcl` (nusystematics' example fcls are not installed by its CMake, hence the source dir) |
| `GENIE_XSEC_TUNE` | `G18_02a_00_000` default; the tune nusystematics' GENIE tools use (`${GENIE_XSEC_TUNE}` in their configs); override per sample |
| `NUISANCE`, `NUISANCE_SRC` | `$NDP_EXTERNAL/nuisance/install` (NUISANCE apps read `$NUISANCE` for `parameters/config.xml` and `data/`), `.../nuisance` |
| PATH / LD_LIBRARY_PATH | `$GENIE_REWEIGHT/bin`, `$NUSYST/bin`, `$NUISANCE/bin` added; libraries `$GENIE_REWEIGHT/lib:$NUWRO/bin:$NUSYST/lib:$NUISANCE/lib` **before** `$ACHILLES/lib` (see 1.7) |

The in-repo GENIE env snapshot `external/genie_env.json` (`pixi run snapshot-genie-env`) only carries
GENIE-related variables; it does not include the NUISANCE/nusystematics ones.

## 1.2 Verify an install

```bash
cd <ndp-platform>
bash .claude/skills/nuisance/scripts/nuisance_check.sh            # installs, programs, libs, env conventions
bash .claude/skills/nuisance/scripts/nuisance_check.sh --smoke    # + pixi run test-nuisance (~10 s)
```

The script re-executes itself under `pixi run` when `$NUISANCE` is not set, so it works from any shell.
`pixi run test-nuisance` (`scripts/test_nuisance.sh`) converts the existing smoke samples and dumps
weights, writing to `external/nuisance-smoke/`; it needs `pixi run test-genie` and `pixi run
test-generators` to have produced `external/genie-smoke/smoke.ghep.root` (50 numu CC on C12 at 3 GeV,
G18_02a_00_000) and `external/generator-smoke/nuwro/smoke.root` (200 events, 3 GeV, C12).

| check (2026-09-13) | expected |
|---|---|
| `PrepareGENIE -i smoke.ghep.root -t "1000060120[1]" -m 3.0` → `nuisflat -i GENIE:…` | `FlatTree_VARS` with 50 events |
| `PrepareNuWroEvents -o … smoke.root` → `nuisflat -i NuWro:…` | `FlatTree_VARS` with 200 events |
| `GenerateSystProviderConfigNuSyst -c $NUSYST_SRC/fcl/ResIso.ToolConfig.fcl` → `DumpConfiguredTweaksNuSyst … -N 20` | `events` tree with 20 entries, `tweak_metadata` with 36; MaRes/MvRes per resonance, weight exactly 1 at the nominal point |
| whole test | 7 s wall |

Manual equivalents: `nuis config --version --features --generators`; `ldd $NUISANCE/bin/nuisflat |
grep nuhepmc_cpputils` must point into `$NUISANCE/lib`; `PrepareGENIE` with no arguments prints its
usage after loading `$NUISANCE/parameters/config.xml` (proves `$NUISANCE` and the libraries are right).

## 1.3 The recipe and how this build deviates from it

`BuildEventGenerators` (README: "All these commands need to be ran on the FNAL gpvm's") sources UPS
products from CVMFS (`setup root v6_28_12 -q e26:p3915:prof`, lhapdf, log4cpp, boost, tbb, sqlite,
hepmc3, …), builds GENIE R-3_06_00 + Reweight R-1_04_00, GiBUU 2025, NuWro, ACHILLES, NEUT (container),
then `build_nusyst.sh` (nusystematics `tags/v02_00_01`, `cmake ../ && make install`) and
`build_nuisance.sh` (NUISANCE default branch, `cmake -DGENIE_ENABLED=ON -DNuWro_ENABLED=ON
-DNEUT_ENABLED=OFF -DGiBUU_ENABLED=ON -DCMAKE_BUILD_TYPE=Debug -DProb3plusplus_ENABLED=ON
-Dnusystematics_ENABLED=ON -DNuHepMC_ENABLED=ON ../`), and `setup_generators.sh` sources
`nusystematics/build/Linux/bin/setup.nusystematics.sh` and `nuisance/build/Linux/setup.sh`.

On the EAF node UPS `setup` fails (`lsb_release: command not found`, known since 2026-09-04), so the
platform builds everything in its pixi environment. Deviations, all deliberate:

| item | recipe | here | why |
|---|---|---|---|
| toolchain | UPS e26 (gcc 12), ROOT 6.28 | pixi gcc 15.2, ROOT 6.40.04, CMake 4.4.3 | UPS unusable; matches the in-repo GENIE/NuWro |
| GENIE + Reweight | R-3_06_00 + R-1_04_00 | R-3_06_02 + R-1_04_02 | current patch releases of the same series |
| nusystematics | v02_00_01 (depends on `jedori0228/systematicstools develop`, a moving branch) | v02_00_07 (pins `NuSystematics/systematicstools v02_00_03`) | reproducible; NUISANCE requires ≥ 2.0.1 |
| NUISANCE | default branch (main) | main, same | commit recorded above |
| build type | Debug | RelWithDebInfo (NUISANCE's own default) | Debug is only a debugging aid |
| install prefixes | `<src>/build/Linux` | `external/<pkg>/install` | platform convention (like ACHILLES) |
| env | source the generated `setup.*.sh` | `activate.sh` exports the same variables | one activation for the whole platform |
| GiBUU | NUISANCE reads its RootTuple output | same code compiled, but the platform's GiBUU has no RootTuple backend | see the `gibuu` skill, section 1.6 |

## 1.4 How the builds work

Prerequisites already in `pixi.toml`: ROOT, gcc/g++/gfortran 15.2 (NUISANCE declares Fortran),
CMake ≥ 3.28, git, GSL, LHAPDF 6, log4cpp, libxml2, HepMC3 ≥ 3.2.5, sqlite (via ROOT), and the three
added on 2026-09-13 for the fhicl-cpp suite / nusystematics: `libboost-devel >=1.75` (got 1.92),
`eigen >=3.4`, `tbb-devel >=2021`. Network access to GitHub is needed at configure time (CPM fetches).

1. **`scripts/build_genie_reweight.sh`** (`pixi run build-genie-reweight`): clones the tag
   `$GENIE_REWEIGHT_VERSION` into `$GENIE_REWEIGHT`, then plain `make -j$NDP_BUILD_JOBS` (default 16)
   with `GOPT_WITH_CXX_USERDEF_FLAGS=-I$ROOTEGPythia6_ROOT/include` (the same TPythia6 trick as
   `build_genie.sh`). Its `src/make/Make.include` includes the Generator's `Make.config` and
   `Make.include`, so GENIE must be configured and built first. No install step; produces
   `lib/libGRw{Fwk,IO,Clc}.so` (+ versioned `-1.04.02-3.06.02` copies, rootmaps, pcm) and `bin/grwght1p`,
   `grwghtnp`, `genie-reweight-config`.
2. **`scripts/build_nusystematics.sh`** (`pixi run build-nusystematics`): clones tag `$NUSYST_VERSION`
   into `$NUSYST_SRC`, then `cmake -S … -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
   -DCMAKE_INSTALL_PREFIX=$NUSYST -DCMAKE_PREFIX_PATH="$CONDA_PREFIX;$ROOTEGPythia6_ROOT"
   -DCMAKE_POLICY_VERSION_MINIMUM=3.5`, build, install. CPM fetches NuHepMC/CMakeModules (its
   `FindGENIE3.cmake` reads `GENIE`, `GENIE_REWEIGHT`, `genie-config --libs`, `PYTHIA6`, `LHAPDF`,
   `LOG4CPP_*`, `LIBXML2_*`, finds ROOTEGPythia6 through `ROOTEGPythia6_ROOT`), systematicstools, and
   Eigen3 (found in pixi). systematicstools builds the fhicl-cpp suite as ExternalProjects
   (`Findfhiclcppstandalone.cmake`: cetmodules 3.22.01 then cetlib-except, hep-concurrency, cetlib,
   fhicl-cpp at FHICLCPP_SUITE_v4_18_01, each `cmake --preset=default`, installed into the same prefix).
   The script exports `CXXFLAGS="-include cassert"` because fhicl-cpp 4.18.01's `DatabaseSupport.cc`
   uses `assert` without the header (gcc 15 error); the ExternalProjects read CXXFLAGS from the
   environment at their configure step.
3. **`scripts/build_nuisance.sh`** (`pixi run build-nuisance`): clones NUISANCE (default branch, or
   `$NUISANCE_REF`) into `$NUISANCE_SRC`, then `cmake … -DCMAKE_INSTALL_PREFIX=$NUISANCE
   -DCMAKE_PREFIX_PATH="$NUSYST;$CONDA_PREFIX;$ROOTEGPythia6_ROOT" -DCMAKE_POLICY_VERSION_MINIMUM=3.5
   -DCMAKE_DISABLE_FIND_PACKAGE_NuHepMC_CPPUtils=TRUE -DGENIE_ENABLED=ON -DNuWro_ENABLED=ON
   -DNEUT_ENABLED=OFF -Dnusystematics_ENABLED=ON -DNuHepMC_ENABLED=ON -DProb3plusplus_ENABLED=ON
   -DT2KReWeight_ENABLED=OFF -DNIWGLegacy_ENABLED=OFF -DNOvARwgt_ENABLED=OFF`, build, install. GENIE
   is found like above (`find_package(GENIE 3.0.0)` → GENIE3), NuWro through `$NUWRO/src/dis/dis_cc.h`
   and `$NUWRO/bin/event1.so`, nusystematics through `find_package(nusystematics 2.0.1)` on the prefix
   path, HepMC3 from pixi; CPM fetches Prob3plusplus 3.10.4 and NuHepMC cpputils (main). Explicitly
   enabled dependencies are *required* by NUISANCE's CMake (configure fails if not found), the three
   `_ENABLED=OFF` ones are skipped instead of searched.

All three scripts are idempotent by existence (clone skipped when `.git` exists; CMake/make rerun
incrementally) and log to `configure.log`, `build.log`, `install.log` in the source directory.

## 1.5 Standalone install outside ndp-platform

Same steps with explicit paths. Requires a GENIE 3 Generator build with `genie-config` on PATH and
`GENIE` exported, NuWro with `event1.so` if `NuWro_ENABLED=ON`, a compiler/ROOT/CMake ≥ 3.21 stack with
Boost ≥ 1.75 (filesystem), SQLite3, TBB ≥ 2020, Eigen3 3.4 (or let CPM fetch it), HepMC3 ≥ 3. Fill every
`xxx`:

```bash
PREFIX=xxx                                   # e.g. a directory on /exp/dune/app
export GENIE=xxx                             # GENIE Generator tree (configured + built; bin/genie-config)
export GENIE_REWEIGHT=$PREFIX/Reweight
git clone --depth 1 --branch R-1_04_02 https://github.com/GENIE-MC/Reweight.git $GENIE_REWEIGHT
( cd $GENIE_REWEIGHT && make -j8 )           # needs $GENIE/src/make/Make.config; add GOPT_WITH_CXX_USERDEF_FLAGS for out-of-ROOT TPythia6 headers

git clone --depth 1 --branch v02_00_07 https://github.com/NuSystematics/nusystematics.git $PREFIX/nusystematics
export CXXFLAGS="-include cassert"           # only needed with gcc >= 15 (fhicl-cpp 4.18.01)
cmake -S $PREFIX/nusystematics -B $PREFIX/nusystematics/build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX=$PREFIX/nusystematics/install -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build $PREFIX/nusystematics/build -j8 && cmake --install $PREFIX/nusystematics/build
export nusystematics_ROOT=$PREFIX/nusystematics/install

git clone https://github.com/NUISANCEMC/nuisance.git $PREFIX/nuisance
export NUWRO=xxx                             # or drop -DNuWro_ENABLED=ON
cmake -S $PREFIX/nuisance -B $PREFIX/nuisance/build -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX=$PREFIX/nuisance/install -DCMAKE_PREFIX_PATH=$nusystematics_ROOT \
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGENIE_ENABLED=ON -DNuWro_ENABLED=ON -DNEUT_ENABLED=OFF \
      -Dnusystematics_ENABLED=ON -DNuHepMC_ENABLED=ON -DProb3plusplus_ENABLED=ON \
      -DT2KReWeight_ENABLED=OFF -DNIWGLegacy_ENABLED=OFF -DNOvARwgt_ENABLED=OFF
cmake --build $PREFIX/nuisance/build -j8 && cmake --install $PREFIX/nuisance/build
source $PREFIX/nuisance/install/setup.sh     # exports NUISANCE and PATH/LD_LIBRARY_PATH (NUISANCE's own script)
source $nusystematics_ROOT/bin/setup.nusystematics.sh
export GENIE_XSEC_TUNE=xxx FHICL_FILE_PATH=.:$nusystematics_ROOT/fcl/fcl:$PREFIX/nusystematics/fcl
```

Add `-DCMAKE_DISABLE_FIND_PACKAGE_NuHepMC_CPPUtils=TRUE` if another NuHepMC cpputils install is
visible on PATH/CMAKE_PREFIX_PATH (1.7). Verify with `nuisance_check.sh --nuisance $PREFIX/nuisance/install
--nusyst $PREFIX/nusystematics/install --reweight $GENIE_REWEIGHT` (no `--smoke` outside the platform).

## 1.6 Updating

- **NUISANCE** tracks `main`: `git -C $NUISANCE_SRC pull` then `pixi run build-nuisance` (incremental);
  for a clean rebuild remove `$NUISANCE_SRC/build` and `$NUISANCE` first. Record the new commit in 1.1.
- **nusystematics / GENIE Reweight** are tags: set `NUSYST_VERSION` / `GENIE_REWEIGHT_VERSION` (in
  `activate.sh` or the environment), move the old source/install directories aside (the scripts skip
  the clone when `.git` exists), rebuild in order Reweight → nusystematics → NUISANCE (NUISANCE links
  both), then `pixi run test-nuisance`.
- Rebuilding GENIE itself invalidates all three (Reweight includes GENIE's Make files; the others link
  GENIE libraries).

## 1.7 Gotchas

- **UPS `setup` is dead on the EAF node**: the recipe's `global_vars.sh` cannot be sourced; use pixi.
- **gcc 15 vs fhicl-cpp 4.18.01**: `'assert' was not declared in this scope` in `DatabaseSupport.cc`;
  fixed by `CXXFLAGS=-include cassert` (build script). Only error in the whole suite.
- **TBB, Boost, Eigen** must exist before configuring nusystematics (`find_package(TBB 2020 REQUIRED)`
  from the fhicl suite finder); now in `pixi.toml`.
- **NuHepMC cpputils clash with ACHILLES**: ACHILLES installed cpputils 0.9.7 (`external/achilles/install`);
  `$ACHILLES/bin` on PATH makes CMake's `find_package` see it, and it lacks `NuHepMC/Reader.hxx` that
  NUISANCE main includes → `CMAKE_DISABLE_FIND_PACKAGE_NuHepMC_CPPUtils=TRUE` forces the CPM fetch. Both
  installs ship `libnuhepmc_cpputils.so` with the same unversioned soname → `activate.sh` lists
  `$NUISANCE/lib` before `$ACHILLES/lib`; ACHILLES still loads its own copy through the `DT_RPATH`
  (`$ORIGIN/../lib`) in its binaries. Symptom when wrong: every NUISANCE program dies with `undefined
  symbol: …NuHepMC…GR8…ReadProcessIdDefinitions…`.
- **`GENIE_XSEC_TUNE` unset** → `GenerateSystProviderConfigNuSyst` aborts with `can't resolve TuneName:
  ${GENIE_XSEC_TUNE}` (the FNAL `genie_xsec` UPS product normally sets it). `activate.sh` defaults it to
  `G18_02a_00_000`; export the sample's tune when it differs.
- **fhicl lookup**: `DumpConfiguredTweaksNuSyst -c file.fcl` searches `FHICL_FILE_PATH` (nusystematics
  ≥ v02_00_06), not the bare path → `.` is on the path in `activate.sh`.
- **`PrepareNuWroEvents` / `PrepareGENIE` refuse to overwrite** an existing output (`file … already
  exists`); remove old outputs first (`test_nuisance.sh` does).
- **NUISANCE's default branch is `main`**, not `master` (the recipe's comment mentions `v2r8`, an old tag).
- **`nuis config --features` omits NuHepMC** even though the input handler is linked; check `ldd`.
- **CMake 4**: `CMAKE_POLICY_VERSION_MINIMUM=3.5` is needed for the older `cmake_minimum_required` of
  CPM-fetched dependencies (Prob3plusplus, Eigen); FetchContent_Populate deprecation warnings (CMP0169)
  from CPM are harmless.
- **`pixi run` prints `WARN Using local manifest …`** on this node because `PIXI_PROJECT_MANIFEST` points
  at the site pixi; harmless.
