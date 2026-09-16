#!/usr/bin/env bash
# Clone (if needed) and build GENIE Reweight $GENIE_REWEIGHT_VERSION in place under external/genie/Reweight.
# Its Make.include pulls in the Generator's Make.config/Make.include, so the in-repo GENIE must be
# configured and built first (pixi run build-genie). No install step: NUISANCE and nusystematics find
# the libraries in $GENIE_REWEIGHT/lib and the headers in $GENIE_REWEIGHT/src.
set -euo pipefail
: "${GENIE:?run via pixi (activate.sh sets GENIE)}"; : "${GENIE_REWEIGHT:?}"; : "${ROOTEGPythia6_ROOT:?}"
{ [ -f "$GENIE/src/make/Make.config" ] && [ -x "$GENIE/bin/genie-config" ]; } || { echo "build GENIE first (pixi run build-genie)" >&2; exit 1; }
JOBS="${NDP_BUILD_JOBS:-16}"
if [ ! -d "$GENIE_REWEIGHT/.git" ]; then
  git clone --depth 1 --branch "$GENIE_REWEIGHT_VERSION" https://github.com/GENIE-MC/Reweight.git "$GENIE_REWEIGHT"
fi
cd "$GENIE_REWEIGHT"
echo "GENIE Reweight $(cat VERSION) ($(git describe --tags --always)) against GENIE $(genie-config --version) | $CXX"
# TPythia6.h & co. come from ROOTEGPythia6 (same trick as build_genie.sh).
export GOPT_WITH_CXX_USERDEF_FLAGS="-I$ROOTEGPythia6_ROOT/include ${GOPT_WITH_CXX_USERDEF_FLAGS:-}"
if ! make -j"$JOBS" > build.log 2>&1; then
  echo "parallel build failed; re-running serially to surface the first real error" >&2
  make > build-serial.log 2>&1 || { grep -m5 -i "error" build-serial.log >&2; exit 1; }
fi
for l in GRwFwk GRwIO GRwClc; do
  [ -f "lib/lib$l.so" ] || { echo "missing lib/lib$l.so (see $GENIE_REWEIGHT/build.log)" >&2; exit 1; }
done
echo "GENIE Reweight built: $(ls lib/*.so | wc -l) libraries, $(ls bin 2>/dev/null | wc -l) binaries in $GENIE_REWEIGHT"
