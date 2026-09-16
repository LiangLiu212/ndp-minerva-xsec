#!/usr/bin/env bash
# Clone (if needed) and CMake-build nusystematics $NUSYST_VERSION into external/nusystematics/install.
# GENIE Generator + Reweight are found through $GENIE / $GENIE_REWEIGHT by the NuHepMC CMake modules
# (fetched by CPM); ROOT, Boost and Eigen3 come from pixi; systematicstools and the standalone
# fhicl-cpp suite (cetlib-except, hep-concurrency, cetlib, fhicl-cpp) are fetched by CPM at configure
# time and installed into the same prefix. CMAKE_POLICY_VERSION_MINIMUM lets CMake 4 accept the old
# cmake_minimum_required of those dependencies.
set -euo pipefail
: "${NUSYST:?run via pixi}"; : "${NUSYST_SRC:?}"; : "${GENIE:?}"; : "${GENIE_REWEIGHT:?}"; : "${ROOTEGPythia6_ROOT:?}"
[ -f "$GENIE_REWEIGHT/lib/libGRwFwk.so" ] || { echo "build GENIE Reweight first (pixi run build-genie-reweight)" >&2; exit 1; }
if [ ! -d "$NUSYST_SRC/.git" ]; then
  git clone --branch "$NUSYST_VERSION" --depth 1 https://github.com/NuSystematics/nusystematics.git "$NUSYST_SRC"
fi
B="$NUSYST_SRC/build"
# fhicl-cpp FHICLCPP_SUITE_v4_18_01 (pinned by systematicstools) uses assert() without <cassert>, which
# gcc 15 rejects; the suite is built as ExternalProjects that read CXXFLAGS from the environment.
export CXXFLAGS="-include cassert ${CXXFLAGS:-}"
echo "nusystematics $(git -C "$NUSYST_SRC" describe --tags --always) | GENIE $(genie-config --version) + Reweight $(cat "$GENIE_REWEIGHT/VERSION") | ROOT $(root-config --version) | $CXX"
cmake -S "$NUSYST_SRC" -B "$B" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="$NUSYST" \
      -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" \
      -DCMAKE_PREFIX_PATH="$CONDA_PREFIX;$ROOTEGPythia6_ROOT" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      > "$NUSYST_SRC/configure.log" 2>&1 || { tail -40 "$NUSYST_SRC/configure.log" >&2; exit 1; }
cmake --build "$B" -j"${NDP_BUILD_JOBS:-16}" > "$NUSYST_SRC/build.log" 2>&1 || { grep -m8 -i -E "error" "$NUSYST_SRC/build.log" >&2; exit 1; }
cmake --install "$B" > "$NUSYST_SRC/install.log" 2>&1
[ -f "$NUSYST/lib/cmake/nusystematics/nusystematicsConfig.cmake" ] && echo "nusystematics built: $NUSYST ($(ls "$NUSYST/bin" | wc -l) programs)" \
  || { echo "nusystematics install incomplete (see $NUSYST_SRC/install.log)" >&2; exit 1; }
