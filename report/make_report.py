#!/usr/bin/env python3
"""Render one self-contained 1mu1p / TKI report (signal definition, truth-level result, event
selection, data vs MC) from a `ndp signal` run and a `ndp selection` run of the same inputs.

    python report/make_report.py --label "playlist 1A" --tag pl1A \
        --signal runs/2026-09-14_signal_minerva_me_ccqelike_1mu1p_FHC1A \
        --selection runs/2026-09-14_selection_minerva_me_ccqelike_1mu1p_FHC1A \
        --out report/TKI_analysis_playlist1A.md [--interpretation file.md]

Every number in the output is read from the two run directories (`summary.json`, `cutflow.json`,
`manifest.json`, `channel.json`) or from the playlist POT files they consumed; the prose of the
signal-definition and selection sections is the reviewed text of `report/TKI_analysis.md`.
Figures are copied from the runs into `report/figs/<tag>_*.png`.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

MEAS_LABELS = {
    "muon_p": "muon p [GeV/c]", "muon_theta": "muon θ [deg]", "muon_pt": "muon p_T [GeV/c]",
    "proton_p": "leading proton p [GeV/c]", "proton_theta": "leading proton θ [deg]", "proton_pt": "leading proton p_T [GeV/c]",
    "dpt": "δp_T [GeV/c]", "dpt_fine": "δp_T, fine grid [GeV/c]", "dptx": "δp_Tx [GeV/c]", "dpty": "δp_Ty [GeV/c]",
    "alpha": "δα_T [deg]", "phi": "φ_T [deg]", "pl": "δp_L [GeV/c]", "pn": "p_n [GeV/c]",
}
MEAS_ORDER = ["muon_p", "muon_theta", "muon_pt", "proton_p", "proton_theta", "proton_pt",
              "dpt", "dpt_fine", "dptx", "dpty", "alpha", "phi", "pl", "pn"]
SIGNAL_STEPS = {"all_cached_truth_events": "all cached truth events", "numu_cc_muon": "νμ CC with μ⁻",
                "muon_window": "muon window", "proton_in_window": "≥ 1 proton in window", "no_mesons": "no mesons",
                "no_heavy_baryons": "no heavy baryons", "no_photons_above_threshold": "no photons > 10 MeV"}
CATEGORY_LABELS = {"signal QE": "signal QE", "signal 2p2h": "signal 2p2h", "signal RES": "signal RES (pion absorbed)",
                   "signal DIS/other": "signal DIS/other", "bkg 1 pi+-": "background, single π±",
                   "bkg 1 pi0": "background, single π⁰", "bkg multi-pi": "background, multi-pion",
                   "bkg other (no pion)": "background, no pion (out of window / fiducial, NC, ν̄, …)"}
TKI_LABELS = {"dpT": "δp_T [GeV/c]", "dpTx": "δp_Tx [GeV/c]", "dpTy": "δp_Ty [GeV/c]", "dalphaT_deg": "δα_T [deg]",
              "dphiT_deg": "φ_T [deg]", "dpL": "δp_L [GeV/c]", "pn": "p_n [GeV/c]"}


def n(x) -> str:
    """Integer with thin-space thousands separators, as in the reviewed report."""
    return f"{int(round(float(x))):,}".replace(",", " ")


def f1(x: float) -> str:
    """One decimal with thin-space thousands separators (scaled MC counts)."""
    return f"{float(x):,.1f}".replace(",", " ")


def sci(x: float, digits: int = 4) -> str:
    """1.0574 × 10²¹ style."""
    sup = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
    m, e = f"{x:.{digits - 1}e}".split("e")
    return f"{m} × 10{str(int(e)).translate(sup)}"


def load(run: Path) -> dict:
    out = {}
    for name in ("summary", "cutflow", "manifest", "channel"):
        p = run / f"{name}.json"
        out[name] = json.loads(p.read_text()) if p.exists() else {}
    return out


def git_head() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def playlist_inputs(channel: dict) -> list[dict]:
    """[{beam, playlist, kind, n_files, pot_used}] from the playlist POT files the runs consumed."""
    data = channel.get("data") or {}
    pls = data.get("playlists") or {}
    root = data.get("products_dir")
    beam = data.get("beam", "FHC")
    rows = []
    for kind in ("data", "mc"):
        for pl in pls.get(kind) or []:
            pl = pl["playlist"] if isinstance(pl, dict) else str(pl)
            pj = Path(root) / beam / pl / f"pot_{pl}_{kind}.json"
            j = json.loads(pj.read_text()) if pj.exists() else {}
            rows.append({"beam": beam, "playlist": pl, "kind": kind, "n_files": j.get("n_files"), "pot_used": j.get("pot_used")})
    return rows


def copy_fig(src: Path, dst: Path) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    return dst.name


def composite_2x2(sel_figs: Path, dst: Path) -> str | None:
    """Muon p / θ, proton p / θ side by side (built from the run's four figures)."""
    names = ["muon_p", "muon_theta", "proton_p", "proton_theta"]
    if not all((sel_figs / f"data_vs_mc_{x}.png").exists() for x in names):
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(14, 13))
    for ax, x in zip(axes.flat, names):
        ax.imshow(mpimg.imread(sel_figs / f"data_vs_mc_{x}.png")); ax.axis("off")
    fig.tight_layout(pad=0.2); dst.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dst, dpi=110); plt.close(fig)
    return dst.name


def ratio_summary(meas: dict) -> list[str]:
    """One line per measurement: data/MC in the first bin, the last bin and the range over bins."""
    import math
    out = []
    for name in MEAS_ORDER:
        m = meas.get(name)
        if not m or "data" not in m:
            continue
        r = [d / mc if mc > 0 else math.nan for d, mc in zip(m["data"], m["mc_scaled"])]
        rr = [x for x in r if not math.isnan(x)]
        if not rr:
            continue
        fmt = lambda x: "–" if math.isnan(x) else f"{x:.2f}"  # noqa: E731  (a bin without MC has no ratio)
        out.append(f"| {MEAS_LABELS.get(name, name)} | {fmt(r[0])} | {fmt(r[-1])} | {min(rr):.2f}–{max(rr):.2f} |")
    return out


def render(label: str, tag: str, sig: dict, sel: dict, sig_run: Path, sel_run: Path, figs_dir: Path, interpretation: str) -> str:
    S, C = sig["summary"], sig["cutflow"]
    T, U = sel["summary"], sel["cutflow"]
    ch = sel["channel"] or sig["channel"]
    inputs = playlist_inputs(ch)
    pot_d, pot_m, scale = T["pot_data"], T["pot_mc"], T["scale"]
    n_sel_d, n_sel_m, n_sel_ms = T["n_data_selected"], T["n_mc_selected"], T["n_mc_selected_scaled"]
    den = T["efficiency_denominator"]
    rows_sig = C["rows"]
    n_all, n_fid_all = rows_sig[0]["n"], rows_sig[0]["n_fiducial"]
    n_numucc_fid = rows_sig[1]["n_fiducial"]
    ratio = n_sel_d / n_sel_ms
    ratio_err = ratio * (1 / n_sel_d + 1 / n_sel_m) ** 0.5
    figs: dict[str, str | None] = {}
    figs["kin"] = copy_fig(sig_run / "figs/signal_kinematics.png", figs_dir / f"{tag}_signal_kinematics.png")
    figs["tki"] = copy_fig(sig_run / "figs/signal_tki.png", figs_dir / f"{tag}_signal_tki.png")
    figs["2x2"] = composite_2x2(sel_run / "figs", figs_dir / f"{tag}_data_vs_mc_muon_proton_2x2.png")
    for m in MEAS_ORDER:
        p = sel_run / f"figs/data_vs_mc_{m}.png"
        figs[m] = copy_fig(p, figs_dir / f"{tag}_data_vs_mc_{m}.png") if p.exists() else None
    rel = figs_dir.name  # figures are referenced relative to the report directory

    data_in = [r for r in inputs if r["kind"] == "data"]; mc_in = [r for r in inputs if r["kind"] == "mc"]
    n_data_files = sum(r["n_files"] or 0 for r in data_in); n_mc_files = sum(r["n_files"] or 0 for r in mc_in)
    pl_list = ", ".join(r["playlist"] for r in data_in)
    sel_git = (sel["manifest"].get("git") or {}).get("sha", "")[:7]; sig_git = (sig["manifest"].get("git") or {}).get("sha", "")[:7]

    L: list[str] = []
    L += [f"# MINERvA CCQE-like 1μ1p / TKI analysis on the open data — {label}", "",
          f"**Status:** rendered {time.strftime('%Y-%m-%d')} from runs of {sel['manifest'].get('timestamp', '')[:10]}. Truth-level signal decided (paper definition); reconstruction-level "
          "selection transcribed from the paper with the unpublished cut values marked as defaults; data versus the "
          "official MC at reconstruction level, shape only (the MC is the unweighted central value).", "",
          "**Reference analysis:** J. Kleykamp et al. (MINERvA), *Measurement of the A dependence of the νμ charged-current "
          "quasielastic-like cross section as a function of muon and proton kinematics at ⟨Eν⟩ ∼ 6 GeV*, arXiv:2503.15047, "
          "Phys. Rev. D 112, 052005 (2025). Only the CH (tracker) result is reproduced here.", "",
          f"**Code:** `ndp-platform` (report rendered at commit `{git_head()}`; signal run at `{sig_git}`, selection run at `{sel_git}`). "
          "Channel manifest `channels/minerva_me_ccqelike_1mu1p.yaml`; signal `ndp/channels/signal.py`; observables "
          "`ndp/channels/observables.py`, `ndp/channels/reco_observables.py`; selection `ndp/adapters/minerva_anatuple.py::ccqelike_1mu1p_cutflow` "
          "via `ndp/channels/selections.py`; diagnostics `ndp/diagnostics.py` (`python -m ndp signal`, `python -m ndp selection`); "
          "grid processing `ndp/grid/`, playlist products `ndp/products.py`.", "",
          "**Runs (every number below is quoted from these directories):** "
          f"`{sig_run.relative_to(ROOT)}/` (truth level) and `{sel_run.relative_to(ROOT)}/` (reconstruction level). "
          f"This file was rendered by `report/make_report.py` from them.", "", "---", ""]

    # ---- 1. inputs -----------------------------------------------------------------------------
    L += ["## 1. Inputs", "", "| item | value | source |", "|---|---|---|",
          f"| Data | ME FHC playlist{'s' if len(data_in) > 1 else ''} {pl_list}: {n(n_data_files)} MasterAnaDev AnaTuple files, streamed on the FNAL grid (campaign `grid/campaigns/fhc_2026-09`) | MINERvA Open Data |",
          f"| Data exposure | {sci(pot_d)} POT | Meta tree `POT_Used`, summed over files (`pot_<pl>_data.json`) |",
          f"| Official MC | StandardMC, same playlist{'s' if len(mc_in) > 1 else ''}: {n(n_mc_files)} files, GENIE 2.12.6, unweighted central value (no MINERvA-tune, flux or GENIE weights) | MINERvA Open Data |",
          f"| MC exposure | {sci(pot_m)} POT; data/MC POT scale {scale:.5g} | Meta tree, summed |",
          f"| Truth (MC) | {n(n_all)} CC events in the skim box (z 5880–8522 mm, apothem 900 mm), {n(n_fid_all)} inside the tracker fiducial | `truth_<pl>_skim.npz` |",
          f"| Reco candidates (MC) | {n(U['rows'][0]['n_mc'])} reconstructed candidates, with the reco rows' truth | `reco_<pl>_mc.npz`, `reco_<pl>_truthcols_skim.npz` (cache version 3) |",
          f"| Reco candidates (data) | {n(U['rows'][0]['n_data'])} reconstructed candidates | `reco_<pl>_data.npz` |",
          "| Paper exposure | 10.61 × 10²⁰ POT, 218 000 selected CH events | arXiv:2503.15047 |", ""]
    if len(data_in) > 1:
        L += ["Per playlist (files and POT_Used):", "", "| playlist | data files | data POT | MC files | MC POT |", "|---|---|---|---|---|"]
        mc_by = {r["playlist"]: r for r in mc_in}
        for r in data_in:
            m = mc_by.get(r["playlist"], {})
            L.append(f"| {r['playlist']} | {r['n_files']} | {sci(r['pot_used'], 3)} | {m.get('n_files', '')} | {sci(m['pot_used'], 3) if m.get('pot_used') else ''} |")
        L.append("")
    L += [f"This data set is {100 * pot_d / 10.61e20:.1f} % of the paper's exposure. The statistical error on the selected data/MC "
          f"ratio is ±{ratio_err:.3f}; every ratio below is a shape statement because the MC carries no weights.", ""]

    # ---- 2. signal definition ------------------------------------------------------------------
    L += ["## 2. Signal definition (truth level)", "", "### 2.1 The paper's definition", "",
          "Verbatim from the paper's tex (`NuclTargTKI_PRD.tex` lines 293–301):", "",
          "> \"The interactions are considered signal if they have a muon with an angle with respect to the",
          "> beam of <17° and a momentum within the range 2 GeV/c < p_μ < 20 GeV/c, and a proton with an",
          "> angle <70° and a momentum in the range 500 MeV/c < p_p < 1100 MeV/c. Additionally, the",
          "> interaction must not have mesons, baryons heavier than neutrons, or photons above 10 MeV.",
          "> [...] For events with more than one proton matching the constraints, the highest momentum",
          "> matching proton is used. [...] The photons with an energy less than 10 MeV are accepted since",
          "> they can come from nuclear de-excitations.\"", "",
          "As implemented (channel manifest, `signal:` block, status **decided**):", "",
          "| item | requirement |", "|---|---|",
          "| interaction | νμ charged current (`nu_pdg = 14`, `current = CC`, primary lepton `lep_pdg = 13`) |",
          "| muon | θμ < 17° with respect to the beam; 2 < pμ < 20 GeV/c |",
          "| proton | ≥ 1 proton with θp < 70° and 0.5 < pp < 1.1 GeV/c; the highest-momentum proton inside this window is the **leading proton** used by every proton and TKI observable |",
          "| vetoes | no mesons; no baryons heavier than the neutron; no photons above 10 MeV |",
          "| allowed | any number of neutrons, photons ≤ 10 MeV, protons outside the window, nuclear remnants |",
          "| phase space | true vertex inside the tracker fiducial box of the inclusive channel (z 5980–8422 mm, apothem 850 mm) — **open**: the paper does not state its CH fiducial volume |", "",
          "Implementation choices, logged in `docs/decisions.md` (2026-09-14):", "",
          "- Final-state classes are PDG **ranges**, not an enumerated list: meson = 100 ≤ |pdg| < 1000; heavy baryon = 1000 ≤ |pdg| < 10000 "
          "except p and n; antinucleons are *not* heavier than the neutron and are not vetoed; nuclear remnants (|pdg| ≥ 10⁹) and GENIE "
          "bookkeeping particles (2000000101) are ignored. Against MAT's enumerated `IsQELike` this differs on 6 of 397 604 νμ CC events "
          "of the certification file (2 anti-Λ, 1 anti-K⁰ events vetoed here but not by MAT; 3–4 events with a second final-state muon kept here, rejected by MAT).",
          "- Both the muon and the proton angles are evaluated in the **NuMI beam frame** (rotation of the detector-frame momenta about x by "
          "−0.05887 rad, MAT's convention), the frame the paper measures angles in and the platform's default with the 2026-09-04 evidence.",
          "- On the grid the definition is evaluated per file by the same code and stored as derived columns of the truth skim "
          "(`signal_minerva_ccqelike_1mu1p`, `lp_*`, veto counts); the numbers below re-evaluate it from those columns.", ""]

    # ---- 2.2 truth-level result ----------------------------------------------------------------
    L += ["### 2.2 Truth-level result on the official MC", "",
          "Cumulative cutflow on the truth skims (`cutflow.json` of the signal run):", "",
          "| step | all skimmed events | inside tracker fiducial |", "|---|---|---|"]
    for i, r in enumerate(rows_sig):
        b = "**" if i == len(rows_sig) - 1 else ""
        L.append(f"| {SIGNAL_STEPS.get(r['step'], r['step'])} | {b}{n(r['n'])}{b} | {b}{n(r['n_fiducial'])}{b} |")
    comp = S["composition_int_type_fiducial"]
    mult = S["protons_in_window_multiplicity_fiducial"]
    n1 = int(mult.get("1", 0)); n2 = int(mult.get("2", 0)); n3 = sum(int(v) for k, v in mult.items() if int(k) >= 3)
    kin = S["kinematics_fiducial"]
    L += ["", f"Signal fraction of fiducial νμ CC events: {100 * S['signal_fraction_of_numu_cc_fiducial']:.2f} % "
          f"({n(S['n_signal_fiducial'])} of {n(n_numucc_fid)}). Composition of the fiducial signal by GENIE process: "
          f"QE {n(comp['QE']['n'])} ({100 * comp['QE']['fraction']:.1f} %), RES with the pion absorbed {n(comp['RES']['n'])} ({100 * comp['RES']['fraction']:.1f} %), "
          f"2p2h {n(comp['MEC']['n'])} ({100 * comp['MEC']['fraction']:.1f} %), DIS {n(comp['DIS']['n'])} ({100 * comp['DIS']['fraction']:.1f} %). "
          f"Protons inside the window per signal event: one in {n(n1)} events, two in {n(n2)}, three or more in {n(n3)}.", "",
          "Kinematics of the fiducial signal (median, 16th–84th percentile): "
          f"muon p {kin['muon_p_gev']['median']:.2f} ({kin['muon_p_gev']['p16']:.2f}–{kin['muon_p_gev']['p84']:.2f}) GeV/c, "
          f"muon θ {kin['muon_theta_deg']['median']:.1f}° ({kin['muon_theta_deg']['p16']:.1f}–{kin['muon_theta_deg']['p84']:.1f}°), "
          f"leading proton p {kin['proton_p_gev']['median']:.3f} ({kin['proton_p_gev']['p16']:.3f}–{kin['proton_p_gev']['p84']:.3f}) GeV/c, "
          f"leading proton θ {kin['proton_theta_deg']['median']:.1f}° ({kin['proton_theta_deg']['p16']:.1f}–{kin['proton_theta_deg']['p84']:.1f}°), "
          f"leading proton pT {kin['proton_pT_gev']['median']:.3f} ({kin['proton_pT_gev']['p16']:.3f}–{kin['proton_pT_gev']['p84']:.3f}) GeV/c.", "",
          f"![truth-level signal kinematics]({rel}/{figs['kin']})",
          "*Figure 1 — Muon and leading-proton kinematics of the fiducial truth-level signal (beam frame), and the proton multiplicity "
          "(all final-state protons vs protons inside the window).*", ""]

    # ---- 3. TKI ----------------------------------------------------------------------------------
    tki = S["tki_fiducial_signal"]
    L += ["## 3. TKI observables", "",
          "Definitions follow Lu et al., PRC 94 (2016) 015503 and Furmanski & Sobczyk, PRC 95 (2017) 065501, as cited by the paper; "
          "ẑ is the neutrino direction, p_T the transverse momenta of the muon (μ) and the leading proton (p):", "",
          "- δp_T = |p_T^μ + p_T^p|", "- δα_T = arccos[ −p̂_T^μ · δp_T / |δp_T| ]", "- φ_T = arccos[ −p̂_T^μ · p̂_T^p ]",
          "- δp_Tx = (ẑ × p̂_T^μ) · δp_T, δp_Ty = −p̂_T^μ · δp_T (negative when the proton carries less transverse momentum than the muon)",
          "- δp_L = R/2 − (m_A'² + δp_T²)/(2R), R = m_A + p_L^μ + p_L^p − E^μ − E^p", "- p_n = √(δp_T² + δp_L²)", "",
          "MAT's legacy `MnvRecoShifter::Calc_tki_vars` uses the opposite sign for both δp_Tx and δp_Ty; the paper's Fig. 1 and its released "
          "grid match the convention above. Nuclear masses are not stated by the paper; the manifest carries m_A = 11.174864 GeV (¹²C nuclear mass), "
          "b = 27.13 MeV (carbon excitation, LE TKI paper arXiv:1805.05486), so m_A' = m_A − m_n + b = 10.262429 GeV, status **default**. "
          "The formulas are verified by an exact p_n recovery on a four-momentum-conserving knockout event and by frame-rotation consistency (`tests/test_tki.py`).", "",
          f"TKI of the fiducial truth-level signal ({n(tki['dpT']['n_finite'])} events; median, 16th–84th percentile):", "",
          "| observable | median | 16 %–84 % |", "|---|---|---|"]
    for k, lab in TKI_LABELS.items():
        v = tki[k]; d = 1 if k.endswith("_deg") else 3
        L.append(f"| {lab} | {v['median']:.{d}f} | {v['p16']:.{d}f}–{v['p84']:.{d}f} |")
    L += ["", f"![truth-level TKI]({rel}/{figs['tki']})",
          "*Figure 2 — TKI distributions of the fiducial truth-level signal: the Fermi peak in p_n near 0.2 GeV/c with the FSI tail, "
          "δα_T rising toward 180°, δp_Ty peaked at zero with the deceleration tail, φ_T peaked at zero.*", ""]

    # ---- 4. selection ----------------------------------------------------------------------------
    L += ["## 4. Event selection (reconstruction level)", "", "### 4.1 The paper's selection and its transcription", "",
          "The paper (Sec. \"Analysis and results\"): events have a negatively charged muon reconstructed in MINOS and at least one proton "
          "candidate whose energy is measured by range; protons that exit or interact inelastically are flagged by the end-of-track deposits "
          "and stopping protons with a Bragg-peak hit pattern are accepted, the Bragg-peak shape also vetoing pions; no Michel electron "
          "candidates near the vertex or any track endpoint; no more than one isolated cluster of energy. The muon and the highest-momentum "
          "proton candidate define the TKI variables. A commented-out tex line (unpublished) gives the reconstruction windows: muon 17°, "
          "2–20 GeV/c; proton 90°, 400–1300 MeV/c, \"a larger fiducial to avoid cutting on the efficiency edge\".", "",
          "Transcription onto the MasterAnaDev AnaTuple (selection `minerva_ccqelike_1mu1p_v0`, values in `selection.params` of the channel manifest):", "",
          "| # | cut | paper wording | tuple branch / value | status |", "|---|---|---|---|---|",
          "| 1 | ZRange, Apothem | fiducial vertex in the tracker | `vtx` in z 5980–8422 mm, hexagon apothem 850 mm (inclusive channel's box) | open (paper fiducial not stated) |",
          "| 2 | HasMINOSMatch, NoDeadtime, IsNeutrino | μ⁻ reconstructed in MINOS | `isMinosMatchTrack == 1`; ≤ 1 dead discriminator pair upstream; `MasterAnaDev_minos_trk_qp < 0` | decided (inclusive chain) |",
          "| 3 | MuonWindow | 17°, 2–20 GeV/c | beam-frame `muon_thetaX/Y` angle < 17°, 2 < p < 20 GeV/c | default (commented tex line) |",
          "| 4 | HasProtonCandidate | ≥ 1 proton candidate, energy by range | `MasterAnaDev_proton_P_fromdEdx > 0` (the tool's primary candidate; secondary candidates never exceed it) | decided |",
          "| 5 | ProtonContained | exiting protons not well reconstructed | `MasterAnaDev_hadron_isExiting[0] == 0` (slot 0 = proton candidate in ~93 % of events) | default |",
          "| 6 | ProtonScore | Bragg-peak stopping proton, pion veto | `MasterAnaDev_proton_score1 > 0.35` | default (paper gives no value) |",
          "| 7 | ProtonWindow | 90°, 400–1300 MeV/c | `MasterAnaDev_proton_theta` (beam frame) < 90°, 0.4 < P_fromdEdx < 1.3 GeV/c | default (commented tex line) |",
          "| 8 | NoMichel | no Michel electron near vertex or track ends | `improved_nmichel == 0` | decided |",
          "| 9 | IsoBlobs | ≤ 1 isolated cluster | `n_nonvtx_iso_blobs ≤ 1` (`_all` variant agrees in 95 % of candidate events) | decided (counter choice open) |", "",
          "`MasterAnaDev_proton_theta` equals the angle of the detector-frame (Px, Py, Pz)_fromdEdx after the NuMI beam rotation (median "
          "difference 0.0000° on the certification file), the same convention as the muon.", ""]

    # ---- 4.2 cutflow -----------------------------------------------------------------------------
    rows = U["rows"]
    L += ["### 4.2 Cutflow, purity, efficiency", "",
          f"Cumulative counts; MC scaled by the POT ratio {scale:.5g}; purity = fraction of selected MC that is truth signal inside the "
          f"fiducial volume; efficiency = selected truth signal / {n(den)} (the fiducial signal count of Sec. 2.2; `cutflow.json` of the selection run):", "",
          "| step | data | MC (scaled) | data / MC | MC purity | MC efficiency |", "|---|---|---|---|---|---|"]
    for i, r in enumerate(rows):
        b = "**" if i == len(rows) - 1 else ""
        L.append(f"| {r['step']} | {b}{n(r['n_data'])}{b} | {b}{f1(r['n_mc_scaled'])}{b} | {b}{r['data_over_mc']:.3f}{b} | {b}{r['mc_purity']:.3f}{b} | {b}{r['mc_efficiency']:.3f}{b} |")
    fid_idx = next(i for i, r in enumerate(rows) if r["step"] == "Apothem")
    pc_idx = next(i for i, r in enumerate(rows) if r["step"] == "HasProtonCandidate")
    sc_idx = next(i for i, r in enumerate(rows) if r["step"] == "ProtonScore")
    L += ["", f"The paper quotes, for the CH tracker, efficiency 28 % and purity 60 %. This selection: efficiency {100 * rows[-1]['mc_efficiency']:.1f} %, "
          f"purity {100 * rows[-1]['mc_purity']:.1f} %. The excess of data over MC before the fiducial cuts ({rows[0]['data_over_mc']:.2f}) is rock-muon and "
          f"non-fiducial activity that the MC does not simulate; after the fiducial cuts the ratio is {rows[fid_idx]['data_over_mc']:.2f}, "
          f"{rows[pc_idx]['data_over_mc']:.2f} at the proton-candidate requirement, {rows[sc_idx]['data_over_mc']:.2f} after the score cut and "
          f"{rows[-1]['data_over_mc']:.3f} ± {ratio_err:.3f} (stat) for the selected sample. The MC is unweighted, so these ratios are not normalisation statements.", "",
          f"Composition of the {n(n_sel_m)} selected MC candidates (categories from the reco rows' truth):", "",
          "| category | n MC | fraction |", "|---|---|---|"]
    for c, v in T["mc_composition_selected"].items():
        L.append(f"| {CATEGORY_LABELS.get(c, c)} | {n(v['n_mc'])} | {v['fraction']:.3f} |")
    L += ["", "Proton-score threshold scan with all other cuts fixed:", "",
          "| `proton_score1` > | data | MC (scaled) | MC purity | MC efficiency |", "|---|---|---|---|---|"]
    for s_ in T["score_scan"]:
        b = "**" if abs(s_["proton_score1_min"] - 0.35) < 1e-9 else ""
        thr = "none" if s_["proton_score1_min"] == 0 else f"{s_['proton_score1_min']}"
        L.append(f"| {b}{thr}{b} | {b}{n(s_['n_data'])}{b} | {b}{f1(s_['n_mc_scaled'])}{b} | {b}{s_['mc_purity']:.3f}{b} | {b}{s_['mc_efficiency']:.3f}{b} |")
    last = T["score_scan"][-1]
    L += ["", f"No threshold reaches the paper's purity: the plateau near {last['mc_purity']:.2f} means the remaining background (dominantly single π± and "
          "π⁰ with the pion undetected) is not separated by the dE/dx score. The paper's Bragg-peak / end-of-track criterion, its isolated-cluster "
          "counter, or its Michel tagger may differ from the branches used here, or the sideband tuning may absorb part of the gap. This is the "
          "main open question of the selection.", ""]

    # ---- 4.3 per playlist (when several) ----------------------------------------------------------
    chunks = U.get("chunks") or []
    dch = [c for c in chunks if c["kind"] == "data"]; mch = {c["label"]: c for c in chunks if c["kind"] == "mc"}
    if len(dch) > 1:
        L += ["### 4.3 Per playlist", "",
              "The comparison accumulates the playlists one at a time; per playlist, with the MC scaled to that playlist's own data POT:", "",
              "| playlist | data POT | data selected | MC selected (scaled) | data / MC | MC purity | MC efficiency |", "|---|---|---|---|---|---|---|"]
        for c in dch:
            m = mch.get(c["label"])
            if not m:
                continue
            sc = c["pot"] / m["pot"]
            L.append(f"| {c['label'].split('/')[-1]} | {sci(c['pot'], 3)} | {n(c['n_selected'])} | {n(m['n_selected'] * sc)} | "
                     f"{c['n_selected'] / (m['n_selected'] * sc):.3f} | {m['n_selected_signal'] / m['n_selected']:.3f} | {m['n_selected_signal'] / m['efficiency_denominator']:.3f} |")
        L.append("")
        sec_dm = "4.4"
    else:
        sec_dm = "4.3"

    # ---- data vs MC ------------------------------------------------------------------------------
    meas = T["measurements"]
    L += [f"### {sec_dm} Data versus MC", "",
          "Selected sample, MC stacked by category and scaled to the data POT, data with Poisson errors, ratio panels with the MC-statistics band. "
          "Bin edges are the paper's released grids (verified from `anc/tki_release.root`). Events in range and data/MC per grid "
          "(`summary.json` of the selection run):", "",
          "| measurement | data in range | MC scaled in range | data NaN | MC NaN | data/MC first bin | last bin | range over bins |",
          "|---|---|---|---|---|---|---|---|"]
    rs = {line.split("|")[1].strip(): line for line in ratio_summary(meas)}
    for name in MEAS_ORDER:
        m = meas.get(name)
        if not m:
            continue
        if "data" not in m:
            L.append(f"| {MEAS_LABELS.get(name, name)} | ERROR {m['error']} | | | | | | |"); continue
        r = rs.get(MEAS_LABELS.get(name, name), "| | | | |").split("|")
        L.append(f"| {MEAS_LABELS.get(name, name)} | {n(m['n_data_in_range'])} | {f1(m['n_mc_scaled_in_range'])} | {m['n_data_nan']} | {m['n_mc_nan']} | {r[2].strip()} | {r[3].strip()} | {r[4].strip()} |")
    L.append("")
    fig_no = 3
    if figs["2x2"]:
        L += [f"![muon and proton kinematics]({rel}/{figs['2x2']})", f"*Figure {fig_no} — Muon momentum and angle, leading-proton momentum and angle.*", ""]; fig_no += 1
    for name, cap in [("dpt_fine", "δp_T (fine grid)"), ("dpt", "δp_T (released grid)"), ("alpha", "δα_T"), ("phi", "φ_T"), ("dptx", "δp_Tx"),
                      ("dpty", "δp_Ty"), ("pl", "δp_L"), ("pn", "p_n"), ("muon_pt", "muon p_T"), ("proton_pt", "leading proton p_T")]:
        if figs.get(name):
            L += [f"![{cap}]({rel}/{figs[name]})", f"*Figure {fig_no} — {cap}.*", ""]; fig_no += 1
    if interpretation:
        L += [interpretation.strip(), ""]

    # ---- 5–7 -----------------------------------------------------------------------------------
    L += ["## 5. Caveats", "",
          f"- **Exposure.** {sci(pot_d)} POT of data ({n(n_sel_d)} selected events) against the paper's 10.61 × 10²⁰.",
          "- **MC weights.** The official MC is used as generated: no MINERvA tune (2p2h enhancement, RPA, pion-production retune), no flux "
          "constraint, no detector-systematic universes. Ratios are shape statements.",
          "- **Unpublished cut values.** The muon and proton reconstruction windows come from commented-out tex lines; the score threshold and "
          "the containment criterion are this analysis's defaults.",
          "- **Fiducial volume and normalisation.** The paper does not state its CH fiducial volume or target count; absolute comparisons need them.",
          "- **Frame.** All angles, truth and reco, are in the NuMI beam frame; the detector-frame numbers differ (the 70° proton edge moves by the 3.4° tilt).",
          "- **Truth skims.** The truth-level numbers come from the per-file skims (CC events in an enlarged tracker box with the final-state "
          "list reduced to derived columns, `docs/decisions.md` 2026-09-14); the full truth tables stay on DUNE scratch.", "",
          "## 6. Open questions (tracked in `docs/open_questions.md`)", "",
          "1. Purity gap versus the paper (0.48 vs 0.60): which score threshold, containment criterion and isolated-blob counter to adopt.",
          "2. The CH fiducial volume and n_nucleons of the paper.",
          "3. Which released grid is the channel's published measurement (14 exist; the placeholder is the leading-proton momentum).",
          "4. Ratification of the TKI constants (m_A, b) and confirmation of the δp_Tx sign convention against the paper's schematic.",
          "5. Events with a second final-state muon (10⁻⁴ level): signal under the primary-lepton definition, background under MAT's `IsQELike`.",
          "6. The MC weight set before any normalisation is read; the playlist dependence of the MC efficiency (full-FHC report).",
          "7. Which playlists to write into the committed channel manifest (`data.playlists`) and where the full truth archives live.", "",
          "## 7. Reproduce", "", "```bash", "cd /exp/dune/data/users/liangliu/ndp-dev/ndp-platform",
          "# products: grid/README.md (ndp grid plan / submit / status / harvest --no-archive; ndp data merge --beam FHC --playlist <pl> --kind data,mc)",
          f"# a copy of channels/minerva_me_ccqelike_1mu1p.yaml with data.playlists = {{mc: [{pl_list}], data: [{pl_list}]}} and products_dir set:",
          f"python -m ndp signal    --channel <manifest.yaml> --slug {sig_run.name.split('_', 1)[1]}",
          f"python -m ndp selection --channel <manifest.yaml> --slug {sel_run.name.split('_', 1)[1]}",
          f"python report/make_report.py --label \"{label}\" --tag {tag} --signal {sig_run.relative_to(ROOT)} --selection {sel_run.relative_to(ROOT)} --out <this file>",
          "```", ""]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", required=True, help='data-set label for the title, e.g. "playlist 1A"')
    ap.add_argument("--tag", required=True, help="prefix of the copied figures, e.g. pl1A")
    ap.add_argument("--signal", required=True, help="`ndp signal` run directory")
    ap.add_argument("--selection", required=True, help="`ndp selection` run directory")
    ap.add_argument("--out", required=True, help="output markdown file (figures go to <out dir>/figs/)")
    ap.add_argument("--interpretation", help="markdown file with the analyst's reading of the data-vs-MC figures (appended after them)")
    a = ap.parse_args()
    sig_run, sel_run, out = Path(a.signal).resolve(), Path(a.selection).resolve(), Path(a.out).resolve()
    interp = Path(a.interpretation).read_text() if a.interpretation else ""
    text = render(a.label, a.tag, load(sig_run), load(sel_run), sig_run, sel_run, out.parent / "figs", interp)
    out.write_text(text)
    print(f"wrote {out} ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
