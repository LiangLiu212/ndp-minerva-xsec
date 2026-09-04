#!/usr/bin/env bash
# Clone (if needed) and CMake-build ACHILLES $ACHILLES_VERSION into external/achilles/install.
# HDF5/zlib/HepMC3/fmt/spdlog/yaml-cpp/HighFive come from pixi (CPM uses local packages when they
# match); the remaining CPM dependencies (NuHepMC cpputils, docopt, yaml-fortran) are fetched at
# configure time. ROOT flux files and Sherpa/BSM are off, as in the user's earlier EAF install.
set -euo pipefail
: "${ACHILLES:?run via pixi}"; : "${ACHILLES_SRC:?}"
if [ ! -d "$ACHILLES_SRC/.git" ]; then
  git clone --branch "$ACHILLES_VERSION" --depth 1 https://github.com/AchillesGen/Achilles.git "$ACHILLES_SRC"
fi
B="$ACHILLES_SRC/build"
cmake -S "$ACHILLES_SRC" -B "$B" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="$ACHILLES" \
      -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_Fortran_COMPILER="$FC" \
      -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DCPM_USE_LOCAL_PACKAGES=ON \
      -DACHILLES_ENABLE_GZIP=ON -DACHILLES_ENABLE_ROOT=OFF -DACHILLES_ENABLE_SHERPA=OFF \
      -DACHILLES_ENABLE_PYTHON=OFF -DACHILLES_ENABLE_TESTING=OFF > "$ACHILLES_SRC/configure.log" 2>&1 \
  || { tail -30 "$ACHILLES_SRC/configure.log" >&2; exit 1; }
cmake --build "$B" -j"${NDP_BUILD_JOBS:-16}" > "$ACHILLES_SRC/build.log" 2>&1 || { grep -m5 -i "error" "$ACHILLES_SRC/build.log" >&2; exit 1; }
cmake --install "$B" > "$ACHILLES_SRC/install.log" 2>&1
[ -x "$ACHILLES/bin/achilles" ] && echo "ACHILLES built: $ACHILLES/bin/achilles ($(git -C "$ACHILLES_SRC" describe --tags --always))" || { echo "ACHILLES install incomplete" >&2; exit 1; }
