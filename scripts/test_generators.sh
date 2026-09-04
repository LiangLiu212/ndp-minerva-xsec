#!/usr/bin/env bash
# Smoke tests for the non-GENIE generators built under external/. Each one runs a tiny job and
# checks that the expected event output exists; GENIE has its own `pixi run test-genie`.
set -euo pipefail
: "${NDP_EXTERNAL:?run via pixi}"
W="$NDP_EXTERNAL/generator-smoke"; mkdir -p "$W"; rc=0

echo "== NuWro =="
if command -v nuwro >/dev/null && [ -f "$NUWRO/bin/event1.so" ]; then
  mkdir -p "$W/nuwro" && cd "$W/nuwro"
  nuwro -o smoke.root -p "number_of_events = 200" -p "beam_type = 0" -p "beam_energy = 3000" -p "beam_particle = 14" \
        -p "target_type = 0" -p "nucleus_p = 6" -p "nucleus_n = 6" -p "random_seed = 7" > nuwro.log 2>&1 || { tail -5 nuwro.log; rc=1; }
  python - <<'PY' || rc=1
import uproot
t = uproot.open("smoke.root")["treeout"]; n = t.num_entries
print(f"NuWro: treeout with {n} events" if n == 200 else f"NuWro: unexpected entry count {n}")
raise SystemExit(0 if n == 200 else 1)
PY
else echo "NuWro not built (pixi run build-nuwro)"; rc=1; fi

echo "== GiBUU =="
if [ -x "$GIBUU/objects/GiBUU.x" ]; then
  mkdir -p "$W/gibuu" && cd "$W/gibuu" && rm -f FinalEvents.dat
  sed "s|__GIBUU_INPUT__|$GIBUU_INPUT|" "$PIXI_PROJECT_ROOT/resources/gibuu/smoke_numu_C12.job" > smoke.job
  "$GIBUU/objects/GiBUU.x" < smoke.job > gibuu.log 2>&1 || { tail -5 gibuu.log; rc=1; }
  if [ -s FinalEvents.dat ]; then echo "GiBUU: FinalEvents.dat with $(grep -vc '^#' FinalEvents.dat) particle rows, $(awk '!/^#/{print $2}' FinalEvents.dat | sort -u | wc -l) events"; else echo "GiBUU: no FinalEvents.dat"; rc=1; fi
else echo "GiBUU not built (pixi run build-gibuu)"; rc=1; fi

echo "== ACHILLES =="
if [ -x "$ACHILLES/bin/achilles" ]; then
  mkdir -p "$W/achilles" && cd "$W/achilles"
  ln -sfn "$ACHILLES_SRC/data" data; cp -f "$ACHILLES_SRC/FormFactors.yml" .
  cat > run.yml <<'YML'
Main:
  NEvents: 200
  HardCuts: false
  EventCuts: false
  DoRotate: false
  RunDecays: False
  Output: {Format: NuHepMC, Name: achilles.hepmc, Zipped: False}
SherpaOptions: !include "data/default/SherpaOptions.yml"
Processes:
  - Leptons: [14, [13]]
Beams:
  - Beam:
      PID: 14
      Beam Params: {Type: Monochromatic, Energy: 1000}
Cascade:
  Run: True
  Interactions: !include "data/default/VirtResInteractions.yml"
  Step: 0.04
  Probability: Cylinder
  InMedium: None
  PotentialProp: False
  Algorithm: Base
NuclearModels:
- NuclearModel:
   Model: QESpectral
   ConfigFile: data/info_C12_pke.data
   SpectralP: data/Spectral_Functions/pke12p_tot.data
   SpectralN: data/Spectral_Functions/pke12n_tot.data
   FormFactorFile: "FormFactors.yml"
   Ward: None
Nuclei:
  - Nucleus: !include "data/default/12C.yml"
Options: !include "data/default/OptionDefaults.yml"
Backend: {Name: Default, Options: []}
YML
  achilles run.yml > achilles.log 2>&1 || true      # achilles exits 0 even on failure: check the log
  if grep -q "Event Run Concluded - Success" achilles.log && [ -s achilles.hepmc ]; then
    echo "ACHILLES: success, $(grep -c '^E ' achilles.hepmc) HepMC3 events"
  else echo "ACHILLES: run did not conclude successfully (see $W/achilles/achilles.log)"; tail -5 achilles.log; rc=1; fi
else echo "ACHILLES not built (pixi run build-achilles)"; rc=1; fi
exit $rc
