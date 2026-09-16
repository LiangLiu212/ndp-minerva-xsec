# Part 1 — Installing GiBUU

Verified on 2026-09-13 on the FNAL EAF node (`jupyter-liangliu-*`, EL9, see the `fnal-eaf-node`
memory) unless another date is given. Sources: the files under `external/gibuu/`, GiBUU's own
`README.md` / `README.Makefile.txt` / `Makefile`, hepforge's `files.txt`, `sha1sum.txt`,
`version.txt`, and the `trac/wiki/compiling` and `trac/wiki/tools` pages.

## 1.1 What this project already has

GiBUU is built inside ndp-platform under `external/` (gitignored, rebuildable from `scripts/`).

| item | value |
|---|---|
| version | Release 2025, patch 5 (April 24, 2026) — `$GIBUU/version.txt` |
| source tree | `external/gibuu/release2025/` (`external/gibuu/release` is a symlink to it) |
| binary | `$GIBUU/objects/GiBUU.x`, 20,498,704 bytes, built 2026-09-04; `$GIBUU/testRun/GiBUU.x` is a symlink to it |
| compiler / flags | pixi gfortran 15.2.0 (conda-forge), `FORT=gfortran`, default `MODE=opt` (-O3), `-fdefault-real-8 -fdefault-double-8 -std=legacy -g -fbacktrace` |
| input tables | `external/gibuu/buuinput` → `buuinput2025/` (what every job card's `path_to_input` must name) |
| tarballs kept | `external/gibuu/release2025.tar.gz` (13,555,950 B), `external/gibuu/buuinput2025.tar.gz` (53,356,096 B) |
| environment | `activate.sh` exports `GIBUU=$NDP_EXTERNAL/gibuu/release2025`, `GIBUU_INPUT=$NDP_EXTERNAL/gibuu/buuinput`, prepends `$GIBUU/objects` to PATH (pixi activation only) |
| build task | `pixi run build-gibuu` → `scripts/build_gibuu.sh` |
| smoke test | `pixi run test-generators` (NuWro, GiBUU, ACHILLES) or `scripts/gibuu_check.sh --smoke` from this skill |
| platform reader | `ndp/adapters/gibuu_finalevents.py` reads `FinalEvents.dat` (`format: gibuu_finalevents` in a model YAML) |
| optional extras | none installed: in-medium tables, RootTuple, HEPMC3event (stub libraries are linked, see 1.4 and 1.6) |

The binary carries `RPATH = <ndp-platform>/.pixi/envs/default/lib` (`readelf -d GiBUU.x`), where
libgfortran.so.5, libquadmath.so.0, libgcc_s.so.1 and libbz2.so.1.0 resolve. It therefore runs
without pixi activation, but stops working if `.pixi/` is removed or the repo is moved. `ldd
GiBUU.x` shows the resolution.

## 1.2 Verify an install

```bash
cd <ndp-platform>
bash .claude/skills/gibuu/scripts/gibuu_check.sh            # tree, libs, tables, hepforge comparison
bash .claude/skills/gibuu/scripts/gibuu_check.sh --smoke    # + numu-C12 smoke job
bash .claude/skills/gibuu/scripts/gibuu_check.sh --gibuu xxx/release2025 --input xxx/buuinput   # another install
```

Healthy result on 2026-09-13: `version.txt` readable; `GiBUU.x` executable; `buuinput/baryon`
present; no `not found` in `ldd`; local tarball sha1s equal to hepforge's; smoke job exit 0.

| smoke job (`resources/gibuu/smoke_numu_C12.job`: numu CC on C12, BNB flux shape, 100 ensembles, 150 time steps) | 2026-09-13 |
|---|---|
| wall time | 25 s |
| `FinalEvents.dat` | 386 events, 1688 particle rows |

The event count is statistical; a value of the same order is fine, zero rows or a non-zero exit
is not. Manual check: the first lines of any GiBUU log print `Version`, `Compiler` and `PATH`.

## 1.3 Download from hepforge

Download page: https://gibuu.hepforge.org/downloads. It is behind an Anubis proof-of-work wall
(version 1.25.0 seen 2026-09-13). Browser-like User-Agents, and Claude's WebFetch tool, get a
"Making sure you're not a bot" page and never the file. Plain `curl` with its default User-Agent
is let through; the direct endpoint is

```
https://gibuu.hepforge.org/downloader?f=<file>        # (downloads?f=<file> 302-redirects here)
```

Do not add `-A "Mozilla/..."`. The same endpoint serves `files.txt` (ls listing), `sha1sum.txt`
and `version.txt`. Files served on 2026-09-13 (all dated Apr 24 on the server):

| file | bytes | sha1 | purpose |
|---|---|---|---|
| release2025.tar.gz | 13,555,950 | 982db881fb54abc9c5b1399545ef3d9322343314 | source; the patch level is only inside (`version.txt`), the name never changes |
| buuinput2025.tar.gz | 53,356,096 | 4d5113eb2c4707a755fafef893820fd81e22eff5 | required input tables |
| buuinput2025_inMed.tar.gz | 405,367,189 | 3d76cddcff6ec4e34acf4e8f18f201a090bd436c | optional in-medium tables (1.6) |
| libraries2025_RootTuple.tar.gz | 23,387 | d0c819e3985ffd80a3968e652267751bd3d1664f | optional ROOT output backend (1.6) |
| libraries2025_HEPMC3.tar.gz | 3,389 | 7d6d50fbc72d29ff300e860b9d9526ee4d360af9 | optional HepMC3 output backend (1.6) |

Fetch and verify:

```bash
cd xxx      # directory that will hold the tarballs
for f in release2025.tar.gz buuinput2025.tar.gz sha1sum.txt version.txt; do
  curl -fsSL --retry 3 -o "$f" "https://gibuu.hepforge.org/downloader?f=$f"
done
sha1sum -c --ignore-missing sha1sum.txt && cat version.txt
```

Both tarballs unpack a versioned directory plus an unversioned symlink (`release → release2025/`,
`buuinput → buuinput2025/`; those symlinks in `external/gibuu/` carry the tarball's own dates).
`scripts/build_gibuu.sh` does the fetch without the sha1 check and reuses a local copy from
`$GIBUU_TARBALLS/<file>` when that variable is set (the 2026-09-04 build reused the nc1p
workspace's tarballs; they are byte-identical to the server's).

Source control: Subversion at `https://gibuu.hepforge.org/repo` (`svn ls
https://gibuu.hepforge.org/svn/` redirects there); not used here. The Trac wiki pages
`trac/wiki/compiling` and `trac/wiki/tools` answer plain curl; `trac/wiki/download` timed out.

## 1.4 Build from source

Prerequisites (GiBUU README + wiki/tools): a Fortran compiler (gfortran ≥ 5.0, ifx/ifort,
nvfortran), GNU make, perl (GiBUU's `Own_Makedepf90.pl` generates the dependencies; `makedepf90`
is the alternative and is not on this node), and libbz2 with headers (the buuinput tables are
bzip2-compressed and decompressed at run time). Here gfortran and bzip2 come from the pixi
environment (`pixi.toml`) and perl 5.32 from `/usr/bin/perl`.

The Makefile takes the first compiler it finds in PATH in the order ifx, ifort, gfortran,
sunf95, pgf95, so always pass `FORT=gfortran` explicitly. Options (from the `Makefile` header):

| variable | values | note |
|---|---|---|
| `FORT` | `ifx`, `ifort`, `gfortran`, `/path/to/compiler` | |
| `MODE` | `opt` (default, -O3), `opt0`…`opt3`, `opt4` (-O3 + vectorize + `-march=native`, gfortran only), `opt5` (auto-parallel, 4 threads), `lto`, `prof`, `callGraph` | platform build uses the default |
| `STATIC` | `0` (default) / `1` | |
| `FPE` | `0`…`3`; `3` = no floating-point-exception trapping (default) | wiki: trapping makes PYTHIA crash |
| `ARGS` | extra compiler flags, e.g. `ARGS="-march=native"` | |
| `withROOT=1`, `withHEPMC3=1` | link the real RootTuple / HEPMC3event libraries | requires 1.6 first |

The command `scripts/build_gibuu.sh` runs:

```bash
cd $GIBUU
make FORT=gfortran -j8 > build.log 2>&1 || make FORT=gfortran > build-serial.log 2>&1
```

`-j` is honoured for the object tree; the serial fallback guards against parallel dependency
generation tripping (not observed here). The 2026-09-04 build compiled 614 source files in
1 min 40 s wall with `-j8` (tarball written 22:47:15, `GiBUU.x` 22:48:55), and `build.log` ends
with `SUCCESS: GiBUU.x generated.` with no lines matching warning/error. The executable lands in
`objects/GiBUU.x`; `testRun/GiBUU.x` is a symlink to it (that is the "in testRun/" of the README).

Stub libraries: the link line always says `-lPDF -lbz2 -lRootTuple -lHEPMC3event`. Without the
extras, `objects/LIB/lib/libRootTuple.a`, `libHEPMC3event.a` and `libPDF.a` are symlinks to stub
archives the build creates itself (`libRootTuplestub.a`, `libHEPMC3eventstub.a`,
`libPDFstub.orig.6225.a`). That is the normal state, not a broken build.

Clean / rebuild (README.Makefile.txt): `make veryclean` removes objects and the executable,
`make superclean` also removes `*.dat`, `fort.*` and editor backups, `make renew` regenerates the
per-directory Makefiles when make complains about missing targets.

## 1.5 Standalone install outside ndp-platform

Same recipe with explicit paths. Fill in every `xxx`:

```bash
PREFIX=xxx                                     # where release2025/ and buuinput2025/ will live
mkdir -p "$PREFIX" && cd "$PREFIX"
for f in release2025.tar.gz buuinput2025.tar.gz sha1sum.txt version.txt; do
  curl -fsSL --retry 3 -o "$f" "https://gibuu.hepforge.org/downloader?f=$f"
done
sha1sum -c --ignore-missing sha1sum.txt
tar xzf release2025.tar.gz && tar xzf buuinput2025.tar.gz
[ -e buuinput ] || ln -s buuinput2025 buuinput
export PATH=xxx:$PATH                          # directory holding gfortran, e.g. <ndp-platform>/.pixi/envs/default/bin
cd release2025 && make FORT=gfortran -j8 > build.log 2>&1; tail -1 build.log
ls -la objects/GiBUU.x
```

Then set `path_to_input='$PREFIX/buuinput'` in every job card (the shipped cards in
`testRun/jobCards/` say `~/GiBUU/buuinput`), and verify with
`gibuu_check.sh --gibuu $PREFIX/release2025 --input $PREFIX/buuinput --smoke`. If gfortran came
from a pixi/conda environment the binary gets an RPATH into that environment's `lib/` (1.1);
keep the environment, or set `LD_LIBRARY_PATH`, or build with a system gfortran.

## 1.6 Optional extras (not installed; procedure from the tarballs' own Makefiles and READMEs, not exercised here)

- **In-medium tables** `buuinput2025_inMed.tar.gz` (405 MB): extra tables for the in-medium
  switches used by `testRun/jobCards/005_Neutrino_*_with_inmed_switches.job` (`&InMedium
  master_piDelta_inmed`, `&XsectionRatios_input flagInMedium`, …). The tarball was not
  downloaded, so how it lays out against `buuinput/` is unverified; run `tar tzf | head` before
  unpacking.
- **RootTuple** (ROOT ntuple output): `libraries2025_RootTuple.tar.gz` unpacks to `libraries2025/`
  plus a `libraries → libraries2025` symlink and must sit *next to* `release2025/` (the Makefile
  sets `LIBDIRSRC = $(ROOTDIR)/../libraries`). Then, with ROOT and cmake on PATH (pixi has ROOT
  6.40): `cd release2025 && make FORT=gfortran buildRootTuple` (cmake + make of RootTuple-master,
  symlinks `libRootTuple.100.a` into `objects/LIB/lib`), then `make FORT=gfortran withROOT=1`.
- **HEPMC3event** (HepMC3 output): same layout from `libraries2025_HEPMC3.tar.gz`;
  `make FORT=gfortran buildHEPMC3event` then `make FORT=gfortran withHEPMC3=1`. Needs HepMC3
  (pixi has hepmc3 ≥ 3.2.5).

The platform reads GiBUU through `FinalEvents.dat`, so none of these is needed for the pipeline.

## 1.7 Updating

GiBUU ships patches by replacing `release2025.tar.gz` in place; only `version.txt` and the sha1
change. `gibuu_check.sh` reports the server's `version.txt` next to the local one and compares
tarball sha1s. To update: fetch into a fresh directory and build there, keep the old tree until
the new binary passes the smoke job, then repoint the `release` symlink or `activate.sh`. In the
platform, move `external/gibuu/release2025` and the tarball aside first, because
`build_gibuu.sh` skips fetch, extract and build whenever their outputs already exist; then
`pixi run build-gibuu`, `pixi run test-generators`, and record the new patch line in 1.1.

## 1.8 Gotchas

- Anubis wall on hepforge (1.3): plain `curl` only; WebFetch and browser agents fail silently.
- `/usr/bin/time` does not exist on the EAF node; time runs with bash's `SECONDS` or `time`.
- The Bash tool's working directory resets between calls; anchor every command with `cd`.
- GiBUU writes its outputs (`FinalEvents.dat`, `*.dat`, `main.run`, `GiBUU_database*`) into the
  current directory; run it from a dedicated run directory.
- `build_gibuu.sh` is idempotent by existence, not by version: a newer tarball is not noticed.
- Make prefers `ifx`/`ifort` over `gfortran` when both are on PATH; always pass `FORT=gfortran`.
- The binary is tied to the pixi environment through its RPATH (1.1).
