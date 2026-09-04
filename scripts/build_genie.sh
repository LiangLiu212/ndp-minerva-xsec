#!/usr/bin/env bash
# Clone (if needed), configure and build GENIE Generator in place under external/genie/Generator,
# against the pixi ROOT / LHAPDF6 / log4cpp / libxml2 / GSL and the ROOTEGPythia6 build.
set -euo pipefail
: "${GENIE:?run via pixi (activate.sh sets GENIE)}"
: "${ROOTEGPythia6_ROOT:?}"
JOBS="${NDP_BUILD_JOBS:-16}"
if [ ! -d "$GENIE/.git" ]; then
  git clone --depth 1 --branch "$GENIE_VERSION" https://github.com/GENIE-MC/Generator.git "$GENIE"
fi
cd "$GENIE"
echo "GENIE $(git describe --tags --always) | compiler $CXX | ROOT $(root-config --version) | LHAPDF $(lhapdf-config --version 2>/dev/null || echo ?)"
# conda-forge ROOT ships no Pythia6: TPythia6.h & co. and libEGPythia6 come from ROOTEGPythia6.
# GENIE's dictionary generation only sees ROOT_INCLUDES (not the user-defined flags) and the
# link line needs the library directory, so patch the two ROOT lines of Make.include once
# (same minimal patch the nc1p workspace used; idempotent).
if ! grep -q "ROOTEGPythia6_ROOT" src/make/Make.include; then
  sed -i 's|^ROOT_INCLUDES  = -I$(shell root-config --incdir)$|ROOT_INCLUDES  = -I$(shell root-config --incdir) -I$(ROOTEGPythia6_ROOT)/include|' src/make/Make.include
  sed -i 's|^ROOT_LIBRARIES = $(shell root-config --glibs) \\$|ROOT_LIBRARIES = $(shell root-config --glibs) -L$(ROOTEGPythia6_ROOT)/lib -Wl,-rpath,$(ROOTEGPythia6_ROOT)/lib \\|' src/make/Make.include
  grep -n "ROOTEGPythia6_ROOT" src/make/Make.include || { echo "Make.include patch did not apply" >&2; exit 1; }
fi
if [ ! -f src/make/Make.config ] || [ "${NDP_RECONFIGURE:-0}" = "1" ]; then
  ./configure \
    --with-compiler=gcc --with-optimiz-level=O2 \
    --enable-pythia6 --with-pythia6-lib="$PYTHIA6" \
    --enable-lhapdf6 --with-lhapdf6-inc="$LHAPDF6_INC" --with-lhapdf6-lib="$LHAPDF6_LIB" \
    --with-libxml2-inc="$LIBXML2_INC" --with-libxml2-lib="$LIBXML2_LIB" \
    --with-log4cpp-inc="$LOG4CPP_INC" --with-log4cpp-lib="$LOG4CPP_LIB" \
    --disable-lhapdf5 2>&1 | tee configure.log
fi
# TPythia6.h & co. come from ROOTEGPythia6, not from ROOT: hand its include dir to every
# compilation (Make.include honours GOPT_WITH_CXX_USERDEF_FLAGS from the environment).
export GOPT_WITH_CXX_USERDEF_FLAGS="-I$ROOTEGPythia6_ROOT/include ${GOPT_WITH_CXX_USERDEF_FLAGS:-}"
echo "building with -j$JOBS (falls back to a serial pass on failure)"
if ! make -j"$JOBS" 2>&1 | tee build.log; then
  echo "parallel build failed; re-running serially to surface the first real error" >&2
  make 2>&1 | tee build-serial.log
fi
ls bin | head -20
echo "GENIE built: $(ls bin | wc -l) binaries, $(ls lib/*.so 2>/dev/null | wc -l) libraries"
