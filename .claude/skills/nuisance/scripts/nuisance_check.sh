#!/usr/bin/env bash
# nuisance_check.sh — verify the GENIE Reweight + nusystematics + NUISANCE installs.
#
# Usage: nuisance_check.sh [--nuisance DIR] [--nusyst DIR] [--reweight DIR] [--smoke] [--no-reexec]
#   --nuisance DIR   NUISANCE install prefix   (default: $NUISANCE)
#   --nusyst DIR     nusystematics install prefix (default: $NUSYST)
#   --reweight DIR   GENIE Reweight tree       (default: $GENIE_REWEIGHT)
#   --smoke          also run the platform smoke test (pixi run test-nuisance)
#   --no-reexec      do not re-execute under `pixi run` when the platform environment is missing
# When $NUISANCE is unset the script re-executes itself under the ndp-platform pixi environment
# (activate.sh provides every variable), so it can be run from any shell.
# Exit status: 0 when every check passes.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"      # physical path: also works via the ndp-dev symlink
PLATFORM="$(cd "$HERE/../../../.." && pwd -P)"                # skill/scripts -> skill -> skills -> .claude -> platform

OPT_NUISANCE=""; OPT_NUSYST=""; OPT_RW=""; SMOKE=0; REEXEC=1
for a in "$@"; do case "$a" in --no-reexec) REEXEC=0;; esac; done
if [ -z "${NUISANCE:-}" ] && [ "$REEXEC" = "1" ] && [ -f "$PLATFORM/pixi.toml" ]; then
  exec pixi run --manifest-path "$PLATFORM/pixi.toml" bash "${BASH_SOURCE[0]}" --no-reexec "$@"
fi
while [ $# -gt 0 ]; do
  case "$1" in
    --nuisance) OPT_NUISANCE=$2; shift 2 ;;
    --nusyst) OPT_NUSYST=$2; shift 2 ;;
    --reweight) OPT_RW=$2; shift 2 ;;
    --smoke) SMOKE=1; shift ;;
    --no-reexec) shift ;;
    -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
NUIS="${OPT_NUISANCE:-${NUISANCE:-}}"; NSY="${OPT_NUSYST:-${NUSYST:-}}"; RW="${OPT_RW:-${GENIE_REWEIGHT:-}}"
fail=0
ok()   { printf '  [ok]   %s\n' "$*"; }
bad()  { printf '  [FAIL] %s\n' "$*"; fail=1; }
info() { printf '  [info] %s\n' "$*"; }
notfound() { ldd "$1" 2>&1 | grep -c "not found" || true; }

echo "== GENIE Reweight: ${RW:-<unset>}"
if [ -n "$RW" ] && [ -f "$RW/VERSION" ]; then
  ok "version $(cat "$RW/VERSION") ($(git -C "$RW" describe --tags --always 2>/dev/null || echo no-git))"
  for l in GRwFwk GRwIO GRwClc; do [ -f "$RW/lib/lib$l.so" ] && ok "lib/lib$l.so" || bad "missing $RW/lib/lib$l.so"; done
  [ -x "$RW/bin/grwght1p" ] && ok "bin/grwght1p present, ldd not-found: $(notfound "$RW/bin/grwght1p")" || bad "missing $RW/bin/grwght1p"
else bad "no GENIE Reweight tree (VERSION file) at '${RW:-}'"; fi
command -v genie-config >/dev/null && info "GENIE $(genie-config --version) at ${GENIE:-?}" || bad "genie-config not on PATH"

echo "== nusystematics: ${NSY:-<unset>}"
if [ -n "$NSY" ] && [ -f "$NSY/lib/cmake/nusystematics/nusystematicsConfig.cmake" ]; then
  ok "nusystematicsConfig.cmake (version $(grep -m1 -o 'nusystematics_VERSION=[0-9.]*' "$NSY/bin/setup.nusystematics.sh" 2>/dev/null | cut -d= -f2)), $(ls "$NSY/bin" | wc -l) programs, $(ls "$NSY/lib" | grep -c '\.so$') shared libraries"
  for b in GenerateSystProviderConfigNuSyst DumpConfiguredTweaksNuSyst fhicl-dump; do
    [ -x "$NSY/bin/$b" ] && ok "bin/$b, ldd not-found: $(notfound "$NSY/bin/$b")" || bad "missing $NSY/bin/$b"
  done
  for l in libnusystematics_systproviders.so libsystematicstools_interface.so libfhiclcpp.so libcetlib.so; do
    [ -f "$NSY/lib/$l" ] || bad "missing $NSY/lib/$l"
  done
else bad "no nusystematics install (nusystematicsConfig.cmake) at '${NSY:-}'"; fi
[ -n "${GENIE_XSEC_TUNE:-}" ] && ok "GENIE_XSEC_TUNE=$GENIE_XSEC_TUNE" || bad "GENIE_XSEC_TUNE unset (nusystematics GENIE tools abort without it)"
case ":${FHICL_FILE_PATH:-}:" in *:.:*) ok "FHICL_FILE_PATH contains '.' ($FHICL_FILE_PATH)";; *) bad "FHICL_FILE_PATH lacks '.' (fcl files in the working directory are not found): '${FHICL_FILE_PATH:-}'";; esac

echo "== NUISANCE: ${NUIS:-<unset>}"
if [ -n "$NUIS" ] && [ -x "$NUIS/bin/nuisflat" ]; then
  ok "$(ls "$NUIS/bin" | wc -l) programs; nuis config --version: $("$NUIS/bin/nuis" config --version 2>/dev/null || echo '?')"
  info "features: $("$NUIS/bin/nuis" config --features 2>/dev/null | tr '\n' ' ')"
  info "generators: $("$NUIS/bin/nuis" config --generators 2>/dev/null | tr '\n' ' ')"
  for b in nuisflat nuiscomp PrepareGENIE PrepareNuWroEvents PrepareGiBUU; do
    [ -x "$NUIS/bin/$b" ] || { bad "missing $NUIS/bin/$b"; continue; }
    n=$(notfound "$NUIS/bin/$b"); [ "$n" = "0" ] && ok "bin/$b resolves all shared libraries" || bad "bin/$b: $n shared librar(ies) not found"
  done
  cpp=$(ldd "$NUIS/bin/nuisflat" 2>/dev/null | awk '/nuhepmc_cpputils/{print $3}')
  if [ -z "$cpp" ]; then info "nuisflat does not load libnuhepmc_cpputils (NuHepMC input handler not built)"
  elif [ "$(cd "$(dirname "$cpp")" && pwd -P)" = "$(cd "$NUIS/lib" && pwd -P)" ]; then ok "nuisflat loads NuHepMC cpputils from \$NUISANCE/lib"
  else bad "nuisflat loads NuHepMC cpputils from $cpp (another install shadows \$NUISANCE/lib; check LD_LIBRARY_PATH order)"; fi
  [ -f "$NUIS/parameters/config.xml" ] && ok "parameters/config.xml present (apps read it through \$NUISANCE)" || bad "missing $NUIS/parameters/config.xml"
  [ "${NUISANCE:-}" = "$NUIS" ] || info "\$NUISANCE (${NUISANCE:-unset}) differs from the checked prefix; the apps use \$NUISANCE"
  [ -x "$NUIS/bin/PrepareGiBUU" ] && info "PrepareGiBUU needs GiBUU RootTuple output (not produced by the platform's GiBUU build)"
else bad "no NUISANCE install (bin/nuisflat) at '${NUIS:-}'"; fi

if [ "$SMOKE" = "1" ]; then
  echo "== smoke: pixi run test-nuisance"
  if [ -x "$PLATFORM/scripts/test_nuisance.sh" ] && [ -n "${NDP_EXTERNAL:-}" ]; then
    SECONDS=0
    if ( cd "$PLATFORM" && bash scripts/test_nuisance.sh > "$PLATFORM/external/nuisance-smoke.log" 2>&1 ); then
      ok "test_nuisance.sh passed in ${SECONDS}s ($(grep -c -E 'flat tree with|events dumped' "$PLATFORM/external/nuisance-smoke.log") checks; outputs in external/nuisance-smoke/)"
    else bad "test_nuisance.sh failed in ${SECONDS}s (see $PLATFORM/external/nuisance-smoke.log)"; tail -5 "$PLATFORM/external/nuisance-smoke.log" | sed 's/^/         /'; fi
  else info "smoke test needs the ndp-platform environment (scripts/test_nuisance.sh); skipped"; fi
fi

echo "== result: $([ "$fail" = "0" ] && echo healthy || echo BROKEN)"
exit $fail
