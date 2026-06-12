# Open questions

Unresolved items the pipeline must settle before (or while) the relevant stage is
implemented. One entry per question; move to **Resolved** with the answer + evidence
when closed.

## Open

### 2026-06-12 — MINOS-efficiency weighter: which beam-intensity input, exactly?

Detail worth remembering for the weights stage: the per-batch `POT_Used_batchN`
leaves in the `Meta` tree are **not** what the MINOS-efficiency reweighter uses —
that one reads per-event beam-intensity branches from the reco tree. The `Meta`
tree is purely exposure bookkeeping.

To settle at Stage-2.3 (weighters): identify the exact reco-tree branch(es) behind
`CVUniverse::GetBatchPOT()` (MAT-MINERvA `MINOSEfficiencyReweighter.h` →
`MinosMuonEfficiencyCorrection::Get(...).GetCorrection(pmu, batchPOT, isFHC)`), and
confirm they are filled in the open-data MasterAnaDev tuples.

## Resolved

(none yet)
