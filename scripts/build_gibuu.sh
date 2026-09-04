#!/usr/bin/env bash
# Fetch GiBUU release 2025 + its input tables and build GiBUU.x with the pixi gfortran.
# Tarballs come from gibuu.hepforge.org; a local copy under $GIBUU_TARBALLS (if set) is reused.
set -euo pipefail
: "${GIBUU:?run via pixi}"; : "${GIBUU_INPUT:?}"
TOP="$(dirname "$GIBUU")"; mkdir -p "$TOP"; cd "$TOP"
fetch() {  # fetch <tarball name>
  local f=$1
  if [ -s "$f" ]; then return; fi
  if [ -n "${GIBUU_TARBALLS:-}" ] && [ -s "$GIBUU_TARBALLS/$f" ]; then cp "$GIBUU_TARBALLS/$f" .; return; fi
  curl -fL --retry 3 -o "$f" "https://gibuu.hepforge.org/downloads?f=$f"
}
fetch release2025.tar.gz; fetch buuinput2025.tar.gz
[ -d release2025 ] || tar -xzf release2025.tar.gz
if [ ! -d buuinput ]; then tar -xzf buuinput2025.tar.gz; [ -d buuinput ] || mv buuinput2025 buuinput; fi
cd release2025
echo "GiBUU $(cat version.txt 2>/dev/null | head -1) | FC=$(command -v gfortran) ($(gfortran --version | head -1))"
if [ ! -x objects/GiBUU.x ]; then
  # GiBUU's own Makefile picks the compiler from PATH ($FORT) and MODE (opt3 default); -j is honoured
  # for the object tree. Fall back to a serial pass if the parallel one trips on a dependency.
  make FORT=gfortran -j"${NDP_BUILD_JOBS:-8}" > build.log 2>&1 || make FORT=gfortran > build-serial.log 2>&1
fi
[ -x objects/GiBUU.x ] && echo "GiBUU.x built: $(ls -la objects/GiBUU.x | awk '{print $5}') bytes" || { echo "GiBUU build failed (see $GIBUU/build.log)" >&2; exit 1; }
[ -d "$GIBUU_INPUT/baryon" ] && echo "buuinput present at $GIBUU_INPUT" || { echo "buuinput missing" >&2; exit 1; }
