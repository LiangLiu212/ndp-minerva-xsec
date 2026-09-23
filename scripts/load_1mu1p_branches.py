#!/usr/bin/env python3
"""Load the AnaTuple branches used by the 1mu + leading proton analysis with uproot.

Standalone (numpy + awkward + uproot only). The branch lists mirror
`ndp/adapters/minerva_anatuple.py` (RECO_BRANCHES, RECO_EXTRA_BRANCHES, RECO_ARRAY_ELEMENTS, TRUTH_*).
Values are returned in the tuple's own units (MeV, mm, radians); nothing is rescaled or NaN-ified.

Usage (repo pixi env: uproot + XRootD):
    pixi run python scripts/load_1mu1p_branches.py <file.root | root://...> [--entry-stop N] [--out file.npz]

From Python:
    from load_1mu1p_branches import load_anatuple
    ev = load_anatuple("MasterAnaDev_mc_AnaTuple_run00110040_Playlist.root", entry_stop=10000)
    ev["reco"]["MasterAnaDev_proton_P_fromdEdx"]      # numpy array, one entry per reco candidate
    ev["truth"]["mc_FSPartPDG"]                        # awkward jagged array, one entry per true event (MC)
    ev["pot"]["pot_used"]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# --- MasterAnaDev tree: reco branches (data and MC) ---------------------------------------------
# One row per reconstructed neutrino candidate. The 1mu1p selection (`minerva_ccqelike_1mu1p_v0`,
# ndp/adapters/minerva_anatuple.py::ccqelike_1mu1p_cutflow) applies twelve cuts in the paper's order:
#   ZRange, Apothem, HasMINOSMatch, NoDeadtime, IsNeutrino, MuonWindow, HasProtonCandidate,
#   ProtonContained, ProtonScore, ProtonWindow, NoMichel, IsoBlobs.
# Cut values below are the channel manifest's (channels/minerva_me_ccqelike_1mu1p.yaml); the ones the
# paper does not print (muon/proton windows, score threshold) are platform defaults, not published numbers.
# Units are the tuple's: MeV, mm, radians. Sentinel -9999 (or -1 for an end point) means "no candidate".

#: used directly by the selection cuts
RECO_SELECTION = (
    # --- event location -------------------------------------------------------------------------
    # Reconstructed interaction vertex [x, y, z, t] in mm.
    #   ZRange:  5980 <= z <= 8422 mm keeps the vertex in the scintillator tracker, away from the
    #            nuclear targets upstream and the calorimeters downstream.
    #   Apothem: |x| < 850 mm and inside the hexagon of apothem 850 mm, i.e. away from the detector edge.
    "vtx",
    # --- muon identification --------------------------------------------------------------------
    # Muon track angles projected on the xz / yz planes, relative to the NuMI BEAM axis (not the
    # detector axis). Combined into the 3D angle theta = acos(1/sqrt(1 + tan^2 thetaX + tan^2 thetaY)).
    #   MuonWindow: theta < 17 deg (muons at larger angles miss MINOS). Also the reco muon angle
    #   observable and the frame reference for the reco TKI variables.
    "muon_thetaX", "muon_thetaY",
    # 1 when the muon track leaving MINERvA is matched to a track in the magnetised MINOS near detector;
    # only matched tracks have a charge and a momentum measurement.
    #   HasMINOSMatch: == 1. This is what defines "a muon" in the analysis.
    "isMinosMatchTrack",
    # Number of dead discriminator pairs along the upstream projection of the muon track. Dead readout
    # there could hide a muon entering from upstream (a rock muon faking a fiducial interaction).
    #   NoDeadtime: <= 1.
    "phys_n_dead_discr_pair_upstream_prim_track_proj",
    # Charge over momentum of the MINOS track fit.
    #   IsNeutrino: < 0 selects a negative track, i.e. mu- from nu_mu (rejects mu+ from anti-nu_mu).
    "MasterAnaDev_minos_trk_qp",
    # Muon 4-vector [px, py, pz, E] in MeV, DETECTOR frame (the platform rotates it about x by the beam
    # angle before combining with the beam-frame angles above).
    #   MuonWindow: 2 < |p| < 20 GeV/c (below 2 GeV/c the MINOS match is unreliable). Also the reco
    #   muon momentum observable and the muon side of the reco TKI variables (dpT, dalphaT, dphiT, ...).
    "MasterAnaDev_leptonE",
    # --- proton identification ------------------------------------------------------------------
    # Momentum (MeV/c) of the leading proton candidate from its range / dE/dx profile under the proton
    # hypothesis; -9999 when the tool found no candidate.
    #   HasProtonCandidate: > 0.   ProtonWindow: 0.4 < p < 1.3 GeV/c.
    #   Also the reco proton momentum observable and the proton side of the reco TKI variables.
    "MasterAnaDev_proton_P_fromdEdx",
    # Angle (rad) of that candidate to the BEAM axis.
    #   ProtonWindow: theta < 90 deg. Also the reco proton angle observable.
    "MasterAnaDev_proton_theta",
    # dE/dx particle-ID score of the candidate: a stopping proton shows a Bragg peak at its end, a pion
    # does not; higher is more proton-like.
    #   ProtonScore: > 0.35 (platform default; the paper prints no number).
    "MasterAnaDev_proton_score1",
    # --- exclusivity vetoes ---------------------------------------------------------------------
    # Number of Michel electrons (pi -> mu -> e decay chain of a stopped pion) near the vertex or a
    # track end point; tags charged pions below tracking threshold.
    #   NoMichel: == 0 (the truth signal has no mesons).
    "improved_nmichel",
    # Number of isolated energy clusters away from the vertex and off every track; pi0 photons leave
    # such clusters, but so can a neutron, hence one is tolerated.
    #   IsoBlobs: <= 1.
    "n_nonvtx_iso_blobs",
)

#: cached for observables and diagnostics (no cut is applied on these)
RECO_EXTRA = (
    # Proton candidate 4-vector components (MeV, DETECTOR frame like the muon), total and kinetic
    # energy: the proton side of the reco TKI variables and the proton momentum/angle cross-checks.
    "MasterAnaDev_proton_Px_fromdEdx", "MasterAnaDev_proton_Py_fromdEdx", "MasterAnaDev_proton_Pz_fromdEdx",
    "MasterAnaDev_proton_E_fromdEdx", "MasterAnaDev_proton_T_fromdEdx",
    # Alternative PID scores (second proton score, pion score): score-threshold studies only.
    "MasterAnaDev_proton_score2", "MasterAnaDev_pion_score1",
    # Candidate track start / end z (mm) and pattern-recognition flag: containment and range cross-checks.
    "MasterAnaDev_proton_startPointZ", "MasterAnaDev_proton_endPointZ", "MasterAnaDev_proton_patternRec",
    # Number of secondary proton candidates (the _sz of the sec_protons array): multi-proton diagnostics.
    "MasterAnaDev_sec_protons_P_fromdEdx_sz",
    # MC only: truth PDG of the particle behind the proton candidate (data: -1). Purity by candidate type.
    "proton_prong_PDG",
    # Isolated-blob alternatives (all blobs, their summed energy), prong and hadron-track multiplicities:
    # the n_nonvtx_iso_blobs vs _all open question and the selection diagnostics.
    "n_nonvtx_iso_blobs_all", "nonvtx_iso_blobs_energy", "n_prongs", "MasterAnaDev_hadron_number",
    # MINOS-only muon momentum (MeV/c): cross-check of the leptonE momentum.
    "MasterAnaDev_minos_trk_p",
    # Recoil (hadronic) energies, four estimators (MeV): calorimetric E_nu and the recoil-E studies.
    "MasterAnaDev_recoil_E", "MasterAnaDev_recoil_passivecorrected", "blob_recoil_E", "recoil_E_polylinecorrected",
    # Total visible energy (MeV).
    "MasterAnaDev_visible_E",
    # The analysis tool's own calorimetric E_nu (MeV), W (MeV), Bjorken x, inelasticity y: reco Q2 and E_nu.
    "MasterAnaDev_E", "MasterAnaDev_W", "MasterAnaDev_x", "MasterAnaDev_y",
)

#: fixed-size array branches of which one element is used: (branch, index)
# MasterAnaDev_hadron_isExiting[0]: 1 if the slot-0 hadron track exits the detector. Slot 0 is the proton
# candidate in ~93 % of events (exploration-repo audit).
#   ProtonContained: == 0 (an exiting proton has no range-based momentum and no Bragg peak).
RECO_ARRAY_ELEMENTS = (("MasterAnaDev_hadron_isExiting", 0),)

# --- truth branches (MC only; present in both the MasterAnaDev and the Truth tree) ---------------
# In MasterAnaDev they describe the true event behind each reco candidate (purity, migration, efficiency
# numerator); in the Truth tree they describe every generated event (efficiency denominator, truth signal).
TRUTH_SCALARS = (
    "mc_incoming",      # neutrino PDG (14 = nu_mu); signal requires nu_mu
    "mc_incomingE",     # true neutrino energy (MeV); flux-weighted E_nu distributions
    "mc_current",       # 1 = charged current; signal requires CC
    "mc_intType",       # MINERvA code: 1 QE, 2 RES, 3 DIS, 4 COH, 8 MEC (2p2h); background/signal by interaction type
    "mc_targetZ", "mc_targetA",   # struck nucleus; signal counts CH (C and H) only
    "mc_Q2", "mc_w",    # true Q2 (MeV^2) and W (MeV); model comparisons
    "mc_primaryLepton", # primary lepton PDG (13 = mu-); signal requires exactly this muon
)
TRUTH_VECTORS = (
    "mc_primFSLepton",  # true muon [px, py, pz, E] MeV, DETECTOR frame -> rotated to the beam frame for
                        # the truth muon window (theta < 17 deg, 2 < p < 20 GeV/c) and the TKI variables
    "mc_vtx",           # true vertex [x, y, z, t] mm; truth fiducial volume (same tracker box as reco)
)
# Jagged final-state particle lists, one entry per particle after FSI (MeV): the leading proton
# (highest p with theta < 70 deg, 0.5 < p < 1.1 GeV/c), the meson / heavy-baryon / photon (> 10 MeV)
# vetoes of the signal definition, and the truth proton side of the TKI variables.
TRUTH_FS = ("mc_nFSPart", "mc_FSPartPDG", "mc_FSPartE", "mc_FSPartPx", "mc_FSPartPy", "mc_FSPartPz")

# --- Meta tree ---------------------------------------------------------------------------------
# One entry per input file: protons on target actually used (POT_Used, summed over files -> the data
# exposure and the data/MC POT scale) and delivered (POT_Total); the entry totals cross-check the trees.
META_BRANCHES = ("POT_Used", "POT_Total", "Total_Reco_Entries", "Total_Truth_Entries")

_TAG = re.compile(r"_(data|mc)_AnaTuple_run0*(\d+)")


def is_mc_file(path: str) -> bool | None:
    """True/False from the AnaTuple file name, None if the name does not follow the open-data pattern."""
    m = _TAG.search(Path(str(path)).name)
    return None if m is None else m.group(1) == "mc"


def _present(tree, wanted) -> tuple[list[str], list[str]]:
    keys = set(tree.keys())
    have = [b for b in wanted if b in keys]
    return have, [b for b in wanted if b not in keys]


def read_meta(f) -> dict:
    """POT and entry totals summed over the Meta tree."""
    tree = f["Meta"]
    have, _ = _present(tree, META_BRANCHES)
    a = tree.arrays(have, library="np")
    out = {"pot_used": float(np.sum(a["POT_Used"])), "pot_total": float(np.sum(a["POT_Total"])),
           "n_meta_entries": int(len(a["POT_Used"]))}
    for k in ("Total_Reco_Entries", "Total_Truth_Entries"):
        if k in a:
            out[k.lower()] = int(np.sum(a[k]))
    return out


def read_reco(f, mc: bool, entry_stop: int | None = None) -> tuple[dict, list[str]]:
    """MasterAnaDev reco branches (+ the matching truth branches for MC) as numpy arrays."""
    tree = f["MasterAnaDev"]
    wanted = list(RECO_SELECTION) + list(RECO_EXTRA) + [b for b, _ in RECO_ARRAY_ELEMENTS]
    if mc:
        wanted += list(TRUTH_SCALARS) + list(TRUTH_VECTORS)
    have, missing = _present(tree, wanted)
    a = tree.arrays(have, library="np", entry_stop=entry_stop)
    for branch, idx in RECO_ARRAY_ELEMENTS:                 # keep the one element the analysis uses
        if branch in a:
            arr = a[branch]
            a[f"{branch}[{idx}]"] = arr[:, idx] if arr.ndim == 2 else np.array([row[idx] for row in arr])
    return a, missing


def read_truth(f, entry_stop: int | None = None) -> tuple[dict, list[str]]:
    """Truth tree (every generated event, MC only): scalars/vectors as numpy, mc_FSPart* as awkward jagged."""
    import awkward as ak
    tree = f["Truth"]
    have, missing = _present(tree, list(TRUTH_SCALARS) + list(TRUTH_VECTORS) + list(TRUTH_FS))
    flat = [b for b in have if b not in TRUTH_FS]
    jag = [b for b in have if b in TRUTH_FS]
    out = dict(tree.arrays(flat, library="np", entry_stop=entry_stop))
    if jag:
        parts = tree.arrays(jag, library="ak", entry_stop=entry_stop)
        for b in jag:
            out[b] = parts[b]                               # ak.Array; ak.flatten / ak.num for CSR form
        if "mc_FSPartPDG" in out:
            out["fs_counts"] = ak.num(out["mc_FSPartPDG"], axis=1).to_numpy()
    return out, missing


def load_anatuple(path: str, entry_stop: int | None = None, mc: bool | None = None, timeout: int = 300) -> dict:
    """Open a local path or root:// URL and load every branch of the 1mu1p analysis."""
    import uproot
    if mc is None:
        mc = is_mc_file(path)
    f = uproot.open(str(path), timeout=int(timeout))        # XRootD wants an integer timeout
    if mc is None:                                          # unknown name: look for truth in the tuple
        mc = "mc_incoming" in f["MasterAnaDev"].keys()
    ev = {"path": str(path), "is_mc": bool(mc), "pot": read_meta(f)}
    ev["reco"], ev["missing_reco"] = read_reco(f, mc, entry_stop)
    if mc and "Truth" in f:
        ev["truth"], ev["missing_truth"] = read_truth(f, entry_stop)
    else:
        ev["truth"], ev["missing_truth"] = {}, []
    return ev


def save_npz(ev: dict, out: str) -> None:
    """Flat npz: reco_<branch>, truth_<branch>; jagged FS branches become truth_<branch>_flat + truth_fs_counts."""
    import awkward as ak
    arrays = {}
    for k, v in ev["reco"].items():
        arrays[f"reco_{k}"] = v
    for k, v in ev["truth"].items():
        if isinstance(v, ak.Array):
            arrays[f"truth_{k}_flat"] = ak.flatten(v).to_numpy()
        else:
            arrays[f"truth_{k}"] = v
    arrays["pot_used"] = np.array(ev["pot"]["pot_used"])
    arrays["pot_total"] = np.array(ev["pot"]["pot_total"])
    np.savez(out, **arrays)


def _summary(ev: dict) -> str:
    lines = [f"{ev['path']}", f"  is_mc: {ev['is_mc']}", f"  POT: {ev['pot']}"]
    n = len(next(iter(ev["reco"].values()))) if ev["reco"] else 0
    lines.append(f"  MasterAnaDev: {n} rows, {len(ev['reco'])} arrays"
                 + (f", missing {ev['missing_reco']}" if ev["missing_reco"] else ""))
    for k in ("MasterAnaDev_proton_P_fromdEdx", "MasterAnaDev_leptonE", "vtx", "MasterAnaDev_hadron_isExiting[0]"):
        if k in ev["reco"]:
            a = ev["reco"][k]
            lines.append(f"    {k}: shape {getattr(a, 'shape', None)} dtype {a.dtype}")
    if ev["truth"]:
        nt = len(ev["truth"]["mc_incoming"]) if "mc_incoming" in ev["truth"] else "?"
        lines.append(f"  Truth: {nt} events, {len(ev['truth'])} arrays"
                     + (f", missing {ev['missing_truth']}" if ev["missing_truth"] else ""))
        if "fs_counts" in ev["truth"]:
            c = ev["truth"]["fs_counts"]
            lines.append(f"    mc_FSPart*: {int(c.sum())} particles, mean {c.mean():.2f} per event")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("path", help="AnaTuple file: local path or root://fndcadoor.fnal.gov:1095/... URL")
    ap.add_argument("--entry-stop", type=int, default=None, help="read only the first N entries of each tree")
    ap.add_argument("--mc", action="store_true", help="force MC handling (default: from the file name)")
    ap.add_argument("--data", action="store_true", help="force data handling")
    ap.add_argument("--out", default=None, help="save the arrays to this .npz")
    args = ap.parse_args(argv)
    mc = True if args.mc else (False if args.data else None)
    ev = load_anatuple(args.path, entry_stop=args.entry_stop, mc=mc)
    print(_summary(ev))
    if args.out:
        save_npz(ev, args.out)
        print(f"  saved -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
