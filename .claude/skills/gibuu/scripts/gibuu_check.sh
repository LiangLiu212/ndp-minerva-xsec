#!/usr/bin/env bash
# gibuu_check.sh — verify a GiBUU install and compare it with what gibuu.hepforge.org serves now.
#
# Usage: gibuu_check.sh [--gibuu DIR] [--input DIR] [--tarballs DIR] [--smoke [JOBCARD]] [--no-server]
#   --gibuu DIR     release tree holding version.txt and objects/GiBUU.x
#                   (default: $GIBUU, else <ndp-platform>/external/gibuu/release2025)
#   --input DIR     buuinput directory (default: $GIBUU_INPUT, else <ndp-platform>/external/gibuu/buuinput)
#   --tarballs DIR  where release2025.tar.gz / buuinput2025.tar.gz live (default: parent of --gibuu)
#   --smoke [CARD]  run a job card in a scratch directory and count the events in FinalEvents.dat
#                   (default card: <ndp-platform>/resources/gibuu/smoke_numu_C12.job; path_to_input is substituted)
#   --no-server     skip the hepforge comparison
# Exit status: 0 when the install is healthy (and the smoke job, if requested, produced events);
# a server mismatch is reported but does not fail the check.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"      # physical path: also works via the ndp-dev symlink
PLATFORM="$(cd "$HERE/../../../.." && pwd -P)"                # skill/scripts -> skill -> skills -> .claude -> platform
SERVER="https://gibuu.hepforge.org/downloader?f="

OPT_GIBUU=""; OPT_INPUT=""; OPT_TARBALLS=""; SMOKE=0; CARD=""; SERVER_CHECK=1
while [ $# -gt 0 ]; do
  case "$1" in
    --gibuu) OPT_GIBUU=$2; shift 2 ;;
    --input) OPT_INPUT=$2; shift 2 ;;
    --tarballs) OPT_TARBALLS=$2; shift 2 ;;
    --smoke) SMOKE=1; shift; if [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; then CARD=$1; shift; fi ;;
    --no-server) SERVER_CHECK=0; shift ;;
    -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

GIBUU_DIR="${OPT_GIBUU:-${GIBUU:-$PLATFORM/external/gibuu/release2025}}"
INPUT_DIR="${OPT_INPUT:-${GIBUU_INPUT:-$PLATFORM/external/gibuu/buuinput}}"
TARBALL_DIR="${OPT_TARBALLS:-$(dirname "$GIBUU_DIR")}"
CARD="${CARD:-$PLATFORM/resources/gibuu/smoke_numu_C12.job}"
BIN="$GIBUU_DIR/objects/GiBUU.x"
fail=0
ok()   { printf '  [ok]   %s\n' "$*"; }
bad()  { printf '  [FAIL] %s\n' "$*"; fail=1; }
info() { printf '  [info] %s\n' "$*"; }

echo "== local install: $GIBUU_DIR"
if [ -r "$GIBUU_DIR/version.txt" ]; then
  LOCAL_VERSION=$(head -1 "$GIBUU_DIR/version.txt"); ok "version.txt: $LOCAL_VERSION"
else LOCAL_VERSION=""; bad "no readable $GIBUU_DIR/version.txt"; fi
if [ -x "$BIN" ]; then
  ok "GiBUU.x: $(stat -c '%s bytes, modified %y' "$BIN" | cut -c1-45)"
  missing=$(ldd "$BIN" 2>&1 | grep -c "not found" || true)
  if [ "$missing" = "0" ]; then ok "ldd: all shared libraries resolve"; else bad "ldd: $missing shared librar(ies) not found"; ldd "$BIN" | grep "not found" | sed 's/^/         /'; fi
  rp=$(readelf -d "$BIN" 2>/dev/null | grep -i -E 'rpath|runpath' | sed -E 's/.*\[(.*)\]/\1/' | head -1)
  [ -n "$rp" ] && info "RPATH: $rp"
else bad "no executable $BIN"; fi
if [ -d "$INPUT_DIR/baryon" ]; then ok "buuinput: $INPUT_DIR ($(ls "$INPUT_DIR" | wc -l) entries)"; else bad "buuinput missing or incomplete: $INPUT_DIR (no baryon/ subdir)"; fi

declare -A LOCAL_SHA
for f in release2025.tar.gz buuinput2025.tar.gz; do
  if [ -s "$TARBALL_DIR/$f" ]; then LOCAL_SHA[$f]=$(sha1sum "$TARBALL_DIR/$f" | cut -c1-40); info "$f: $(stat -c %s "$TARBALL_DIR/$f") bytes, sha1 ${LOCAL_SHA[$f]}"
  else info "$f: not kept in $TARBALL_DIR"; fi
done

if [ "$SERVER_CHECK" = "1" ]; then
  echo "== hepforge (plain curl; the site's bot wall blocks browser user agents)"
  SV=$(curl -fsS --max-time 20 "${SERVER}version.txt" 2>/dev/null | head -1 || true)
  if [ -z "$SV" ]; then info "server unreachable or blocked; skipping comparison"
  else
    if [ "$SV" = "$LOCAL_VERSION" ]; then ok "server version.txt matches: $SV"; else info "server version.txt: $SV  (local: ${LOCAL_VERSION:-none}) -> newer/different release available"; fi
    SUMS=$(curl -fsS --max-time 20 "${SERVER}sha1sum.txt" 2>/dev/null || true)
    for f in "${!LOCAL_SHA[@]}"; do
      ssha=$(printf '%s\n' "$SUMS" | awk -v f="$f" '$2==f{print $1}')
      if [ -z "$ssha" ]; then info "$f: not listed in server sha1sum.txt"
      elif [ "$ssha" = "${LOCAL_SHA[$f]}" ]; then ok "$f: sha1 identical to the server's"
      else info "$f: sha1 differs from the server's ($ssha) -> re-fetch to update"; fi
    done
  fi
fi

if [ "$SMOKE" = "1" ]; then
  echo "== smoke job: $CARD"
  if [ ! -x "$BIN" ] || [ ! -r "$CARD" ]; then bad "cannot run smoke job (binary or job card missing)"
  else
    W=$(mktemp -d "${TMPDIR:-/tmp}/gibuu_smoke.XXXXXX")
    sed -e "s|__GIBUU_INPUT__|$INPUT_DIR|" -e "s|path_to_input *= *'[^']*'|path_to_input='$INPUT_DIR'|" "$CARD" > "$W/smoke.job"
    SECONDS=0
    ( cd "$W" && "$BIN" < smoke.job > gibuu.log 2>&1 ); rc=$?
    if [ "$rc" != "0" ]; then bad "GiBUU.x exited $rc after ${SECONDS}s (log: $W/gibuu.log)"; tail -5 "$W/gibuu.log" | sed 's/^/         /'
    elif [ -s "$W/FinalEvents.dat" ]; then
      rows=$(grep -vc '^#' "$W/FinalEvents.dat"); evs=$(awk '!/^#/{print $1" "$2}' "$W/FinalEvents.dat" | sort -u | wc -l)
      if [ "$evs" -gt 0 ]; then ok "exit 0 in ${SECONDS}s: FinalEvents.dat with $evs events, $rows particle rows (dir: $W)"; else bad "FinalEvents.dat has no events (dir: $W)"; fi
    else bad "exit 0 in ${SECONDS}s but no FinalEvents.dat (dir: $W)"; fi
  fi
fi

echo "== result: $([ "$fail" = "0" ] && echo healthy || echo BROKEN)"
exit $fail
