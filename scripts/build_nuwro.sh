#!/usr/bin/env bash
# Clone (if needed) and build NuWro in place under external/nuwro with the pixi ROOT and the
# ROOTEGPythia6 Pythia6 (NuWro's Makefile already honours $ROOTEGPythia6_ROOT).
set -euo pipefail
: "${NUWRO:?run via pixi}"; : "${ROOTEGPythia6_ROOT:?}"
[ -f "$ROOTEGPythia6_ROOT/lib/libEGPythia6.so" ] || { echo "build ROOTEGPythia6 first (pixi run build-pythia6)" >&2; exit 1; }
NUWRO_REF="${NUWRO_REF:-master}"
if [ ! -d "$NUWRO/.git" ]; then git clone https://github.com/NuWro/nuwro.git "$NUWRO"; git -C "$NUWRO" checkout -q "$NUWRO_REF"; fi
cd "$NUWRO"; echo "NuWro $(git describe --tags --always) | ROOT $(root-config --version) | $(g++ --version | head -1)"
make -j"${NDP_BUILD_JOBS:-16}" 2>&1 | tee build.log | grep -E "error|Error" || true
ls bin | tr '\n' ' '; echo
[ -x bin/nuwro ] && [ -f bin/event1.so ] && echo "NuWro built" || { echo "NuWro build incomplete (see $NUWRO/build.log)" >&2; exit 1; }
