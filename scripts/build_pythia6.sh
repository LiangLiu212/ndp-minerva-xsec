#!/usr/bin/env bash
# Build luketpickering/ROOTEGPythia6 (ROOT's EGPythia6 library, removed from ROOT 6.30, plus a
# built-in Pythia6) with the pixi compilers into external/ROOTEGPythia6/install.
set -euo pipefail
: "${NDP_EXTERNAL:?run via pixi (activate.sh sets NDP_EXTERNAL)}"
SRC="$NDP_EXTERNAL/ROOTEGPythia6"
if [ ! -d "$SRC/.git" ]; then
  git clone https://github.com/luketpickering/ROOTEGPythia6.git "$SRC"
fi
if [ -f "$ROOTEGPythia6_ROOT/lib/libEGPythia6.so" ] && [ -f "$ROOTEGPythia6_ROOT/lib/libPythia6.so" ]; then
  echo "ROOTEGPythia6 already installed at $ROOTEGPythia6_ROOT"; exit 0
fi
# -march=native is ROOTEGPythia6's default for Pythia6; keep binaries portable across the
# EAF/grid host families by compiling for a generic x86-64 level instead.
sed -i 's/target_compile_options(Pythia6 PRIVATE -march=native)/target_compile_options(Pythia6 PRIVATE -march=x86-64-v2)/' "$SRC/CMakeLists.txt"
cmake -S "$SRC" -B "$SRC/build" -DROOTEGPythia6_Pythia6_BUILTIN=ON -DCMAKE_INSTALL_PREFIX="$ROOTEGPythia6_ROOT" \
      -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_Fortran_COMPILER="$FC"
cmake --build "$SRC/build" -j"${NDP_BUILD_JOBS:-16}"
cmake --install "$SRC/build"
ls "$ROOTEGPythia6_ROOT/lib"
