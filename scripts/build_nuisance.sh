#!/usr/bin/env bash
# Clone (if needed) and CMake-build NUISANCE (default branch, or $NUISANCE_REF) into external/nuisance/install
# against the in-repo GENIE (+Reweight), NuWro, nusystematics and the pixi ROOT / HepMC3. Generators are
# located through the environment (GENIE, GENIE_REWEIGHT, NUWRO, nusystematics_ROOT); NuHepMC cpputils
# and Prob3plusplus are fetched by CPM (the find_package for cpputils is disabled on purpose: ACHILLES installed
# cpputils 0.9.7 on PATH, which lacks NuHepMC/Reader.hxx that NUISANCE main needs). NEUT is off (not installed). GiBUU input support is compiled in
# unconditionally by NUISANCE (it reads GiBUU's RootTuple ROOT output, which needs GiBUU built withROOT=1).
set -euo pipefail
: "${NUISANCE:?run via pixi}"; : "${NUISANCE_SRC:?}"; : "${GENIE:?}"; : "${GENIE_REWEIGHT:?}"; : "${NUWRO:?}"; : "${NUSYST:?}"; : "${ROOTEGPythia6_ROOT:?}"
[ -f "$NUSYST/lib/cmake/nusystematics/nusystematicsConfig.cmake" ] || { echo "build nusystematics first (pixi run build-nusystematics)" >&2; exit 1; }
[ -f "$NUWRO/bin/event1.so" ] || { echo "build NuWro first (pixi run build-nuwro)" >&2; exit 1; }
if [ ! -d "$NUISANCE_SRC/.git" ]; then
  git clone https://github.com/NUISANCEMC/nuisance.git "$NUISANCE_SRC"        # default branch (main) unless NUISANCE_REF says otherwise
  [ -n "${NUISANCE_REF:-}" ] && git -C "$NUISANCE_SRC" checkout -q "$NUISANCE_REF"
fi
B="$NUISANCE_SRC/build"
echo "NUISANCE $(git -C "$NUISANCE_SRC" describe --tags --always) | GENIE $(genie-config --version) | NuWro $(git -C "$NUWRO" describe --tags --always) | ROOT $(root-config --version) | $CXX"
cmake -S "$NUISANCE_SRC" -B "$B" -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX="$NUISANCE" \
      -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_Fortran_COMPILER="$FC" \
      -DCMAKE_PREFIX_PATH="$NUSYST;$CONDA_PREFIX;$ROOTEGPythia6_ROOT" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DCMAKE_DISABLE_FIND_PACKAGE_NuHepMC_CPPUtils=TRUE \
      -DGENIE_ENABLED=ON -DNuWro_ENABLED=ON -DNEUT_ENABLED=OFF -Dnusystematics_ENABLED=ON \
      -DNuHepMC_ENABLED=ON -DProb3plusplus_ENABLED=ON \
      -DT2KReWeight_ENABLED=OFF -DNIWGLegacy_ENABLED=OFF -DNOvARwgt_ENABLED=OFF \
      > "$NUISANCE_SRC/configure.log" 2>&1 || { tail -40 "$NUISANCE_SRC/configure.log" >&2; exit 1; }
cmake --build "$B" -j"${NDP_BUILD_JOBS:-16}" > "$NUISANCE_SRC/build.log" 2>&1 || { grep -m8 -i -E "error" "$NUISANCE_SRC/build.log" >&2; exit 1; }
cmake --install "$B" > "$NUISANCE_SRC/install.log" 2>&1
[ -x "$NUISANCE/bin/nuisflat" ] && echo "NUISANCE built: $NUISANCE ($(ls "$NUISANCE/bin" | wc -l) programs)" \
  || { echo "NUISANCE install incomplete (see $NUISANCE_SRC/install.log)" >&2; exit 1; }
