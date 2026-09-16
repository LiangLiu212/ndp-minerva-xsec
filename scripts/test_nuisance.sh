#!/usr/bin/env bash
# Smoke test for the NUISANCE / nusystematics stack built under external/: converts the GENIE and NuWro
# smoke samples (from `pixi run test-genie` and `pixi run test-generators`) into NUISANCE flat trees, and
# runs a nusystematics GENIE-reweight dump (ResIso provider: MaRes/MvRes per resonance) on the GENIE
# sample. Outputs go to external/nuisance-smoke/.
set -euo pipefail
: "${NDP_EXTERNAL:?run via pixi}"; : "${NUISANCE:?}"; : "${NUSYST:?}"; : "${NUSYST_SRC:?}"
W="$NDP_EXTERNAL/nuisance-smoke"; mkdir -p "$W"; cd "$W"; rc=0
rm -f ./*.root ./*.fcl ./*.log          # the Prepare* tools refuse to overwrite existing output files
GHEP="$NDP_EXTERNAL/genie-smoke/smoke.ghep.root"          # 50 numu CC on C12 at 3 GeV, G18_02a_00_000
NUWRO_OUT="$NDP_EXTERNAL/generator-smoke/nuwro/smoke.root"  # 200 numu on C12 at 3 GeV
count() {  # count <file> [tree]: entries of <tree>, or of the largest TTree in the file
python - "$1" "${2:-}" <<'PY'
import sys, uproot
u = uproot.open(sys.argv[1]); want = sys.argv[2]
trees = [(k.split(";")[0], u[k].num_entries) for k in u.keys() if hasattr(u[k], "num_entries")]
print(dict(trees).get(want, 0) if want else (max(t[1] for t in trees) if trees else 0))
PY
}

echo "== GENIE -> PrepareGENIE -> nuisflat =="
if [ -f "$GHEP" ]; then
  PrepareGENIE -i "$GHEP" -t "1000060120[1]" -m 3.0 -o genie.prep.root > prep_genie.log 2>&1 || { tail -5 prep_genie.log; rc=1; }
  nuisflat -i GENIE:genie.prep.root -o genie.flat.root > flat_genie.log 2>&1 || { tail -5 flat_genie.log; rc=1; }
  n=$(count genie.flat.root); [ "$n" = "50" ] && echo "GENIE: flat tree with $n events" || { echo "GENIE: unexpected flat-tree entries: $n"; rc=1; }
else echo "GENIE smoke sample missing ($GHEP): run pixi run test-genie"; rc=1; fi

echo "== NuWro -> PrepareNuWroEvents -> nuisflat =="
if [ -f "$NUWRO_OUT" ]; then
  PrepareNuWroEvents -o nuwro.prep.root "$NUWRO_OUT" > prep_nuwro.log 2>&1 || { tail -5 prep_nuwro.log; rc=1; }
  nuisflat -i NuWro:nuwro.prep.root -o nuwro.flat.root > flat_nuwro.log 2>&1 || { tail -5 flat_nuwro.log; rc=1; }
  n=$(count nuwro.flat.root); [ "$n" = "200" ] && echo "NuWro: flat tree with $n events" || { echo "NuWro: unexpected flat-tree entries: $n"; rc=1; }
else echo "NuWro smoke sample missing ($NUWRO_OUT): run pixi run test-generators"; rc=1; fi

echo "== nusystematics: ResIso provider headers + weights on the GENIE sample =="
if [ -f "$GHEP" ]; then
  export GENIE_XSEC_TUNE=G18_02a_00_000          # the tune the smoke sample was generated with
  GenerateSystProviderConfigNuSyst -c "$NUSYST_SRC/fcl/ResIso.ToolConfig.fcl" -o resiso.headers.fcl > gen_resiso.log 2>&1 || { tail -5 gen_resiso.log; rc=1; }
  DumpConfiguredTweaksNuSyst -c resiso.headers.fcl -i "$GHEP" -N 20 -o resiso.tweaks.root > dump_resiso.log 2>&1 || { tail -5 dump_resiso.log; rc=1; }
  n=$(count resiso.tweaks.root events); nom=$(grep -c "variation\[3\] = 0 calculated weight: 1$" dump_resiso.log || true)
  [ "$n" = "20" ] && [ "$nom" -gt 0 ] && echo "nusystematics: $n events dumped, $nom nominal-point weights exactly 1" || { echo "nusystematics: unexpected output (entries=$n, nominal-weight lines=$nom)"; rc=1; }
fi

[ "$rc" = "0" ] && echo "NUISANCE smoke tests passed (outputs in $W)" || echo "NUISANCE smoke tests FAILED (see $W)"
exit $rc
