"""MnvTune v1 CV weights #2-#5: non-res pi, 2p2h, RPA, MINOS efficiency.

Python ports of the MAT-MINERvA weighters (see docs/cv_reweight.md for the
physics; per-function provenance below). Together with xsec.flux.FluxCV
(weight #1) they reproduce Model::GetWeight at CV:

    truth-only product (ALL MC fills, incl. efficiency denominator):
        flux x nonres_pi x 2p2h x rpa
    reco fills additionally:
        x minos_efficiency

All functions are vectorized; kinematic inputs in the units noted.
External files come from the unpacked tarball under data/flux/.
"""
from functools import lru_cache
from pathlib import Path

import numpy as np
import uproot

from .constants import CC_CURRENT, NONRES_PI_WEIGHT, NU_MU_PDG

REPO_ROOT = Path(__file__).resolve().parents[1]
REWEIGHT_DIR = REPO_ROOT / "data" / "flux" / "MParamFiles" / "data" / "Reweight"

# GENIE interaction types (mc_intType) used by the weighters.
INT_TYPE_QE = 1
INT_TYPE_MEC = 8


# ---------------------------------------------------------------------------
# #2 non-resonant pion reduction — GENIEReweighter(true, false)
# ---------------------------------------------------------------------------
def _knob_element(arr, idx=2):
    """Element [idx] of a GENIE knob array per event; tolerates both the
    fixed (n, k) 2D layout and uproot's jagged object-array layout."""
    a = np.asarray(arr)
    if a.dtype == object:
        return np.array([row[idx] for row in a], dtype=np.float64)
    return a[:, idx].astype(np.float64)


def nonres_pi_weight(rvn1pi, rvp1pi):
    """w = 0.43 for GENIE-tagged non-res single-pi events, else 1.

    Inputs are the truth_genie_wgt_Rv{n,p}1pi knob arrays; the tag is
    element [2] < 1 on either (GenieSystematics.cxx:368-373).
    """
    tagged = (_knob_element(rvn1pi) < 1.0) | (_knob_element(rvp1pi) < 1.0)
    return np.where(tagged, NONRES_PI_WEIGHT, 1.0)


# ---------------------------------------------------------------------------
# #3 low-recoil 2p2h enhancement — LowRecoil2p2hReweighter, mode 0 (CV)
# ---------------------------------------------------------------------------
class TwoP2HWeight:
    """w = 1 + correlated 2D Gaussian in true (q0, q3) [GeV] for MEC events.

    CV (variation 0) applies ONLY to mc_intType==8 on nuclei (mc_targetZ>=2);
    QE stays 1 at CV (MnvTuneSystematics.cxx:19-60). Gaussian and the
    "+1" convention: weight_2p2h.{h,cxx} (Gaussian2Dplusone). Fit parameters:
    MParamFiles/data/Reweight/fit-mec-2d-noScaleDown-penalty00300-best-fit.txt
    (order: norm, meanq0, meanq3, sigmaq0, sigmaq3, corr).
    """

    def __init__(self, path=None):
        path = Path(path or REWEIGHT_DIR / "fit-mec-2d-noScaleDown-penalty00300-best-fit.txt")
        vals = [float(x) for x in path.read_text().split()][:6]
        (self.norm, self.meanq0, self.meanq3,
         self.sigmaq0, self.sigmaq3, self.corr) = vals

    def gaussian2d(self, q0, q3):
        z = ((q0 - self.meanq0) ** 2 / self.sigmaq0 ** 2
             + (q3 - self.meanq3) ** 2 / self.sigmaq3 ** 2
             - 2 * self.corr * (q0 - self.meanq0) * (q3 - self.meanq3)
             / (self.sigmaq0 * self.sigmaq3))
        return self.norm * np.exp(-0.5 * z / (1 - self.corr ** 2))

    def weight(self, q0_gev, q3_gev, mc_intType, mc_targetZ):
        q0 = np.asarray(q0_gev, dtype=np.float64)
        q3 = np.asarray(q3_gev, dtype=np.float64)
        applies = ((np.asarray(mc_intType) == INT_TYPE_MEC)
                   & (np.asarray(mc_targetZ) >= 2))
        return np.where(applies, 1.0 + self.gaussian2d(q0, q3), 1.0)


# ---------------------------------------------------------------------------
# #5 RPA suppression — RPAReweighter (CV), weightRPA::getWeightInternal
# ---------------------------------------------------------------------------
# q0 offsets for useNX=true (GENIE 2.12.X era), nu mode, rpaMat=true
# (weightRPA.cxx:30-53); default 27 for nuclei without a dedicated table.
RPA_Q0_OFFSET_NU = {6: 27, 8: 25, 26: 39, 82: 35}
RPA_Q0_OFFSET_DEFAULT = 27
# Z -> ratio file (weightRPA.cxx:367-418); other Z>=6 fall back to 12C.
RPA_FILE_NU = {6: "outNievesRPAratio-nu12C-20GeV-20170202.root",
               8: "outNievesRPAratio-nu16O-20GeV-20170202.root",
               26: "outNievesRPAratio-nu56Fe-20GeV-20170202.root",
               82: "outNievesRPAratio-nu208Pb-20GeV-20170202.root"}


class RPAWeight:
    """Valencia RPA / GENIE ratio at true (q0, q3) [GeV] for true-QE on Z>=6.

    Port of weightRPA::getWeightInternal (weightRPA.cxx:55-129), CV only:
    MeV-bin index lookup into hrelratio (3000x3000 over [0,3]^2 GeV, axes
    x=q3, y=q0) with the documented traps: >=3 GeV clamp to bin 2999;
    q0<0.018 -> row 18; w<=0.001 -> 1; (q0<0.15 & w>0.9) -> re-look at
    q3bin+150; Q2=q3^2-q0^2 in (3,9) -> hQ2relratio at FindBin(Q2), >=9 -> 1;
    final sanity w outside [0.001, 2] -> 1. Gating: mc_intType==1 and
    mc_targetZ>=6 (MnvTuneSystematics.cxx:103-105).
    """

    def __init__(self, nu_pdg=NU_MU_PDG):
        if nu_pdg != NU_MU_PDG:
            raise NotImplementedError("only nu (FHC) tables wired")
        self._tables = {}

    @lru_cache(maxsize=8)
    def _load(self, fname):
        with uproot.open(REWEIGHT_DIR / fname) as f:
            rel = f["hrelratio"].values(flow=True)          # (3002, 3002): [q3bin, q0bin]
            q2vals = f["hQ2relratio"].values(flow=True)     # (10002,)
            q2edges = f["hQ2relratio"].axis(0).edges()
        return rel, q2vals, q2edges

    def _weight_one_z(self, q0, q3, z):
        fname = RPA_FILE_NU.get(int(z), RPA_FILE_NU[6])
        offset = RPA_Q0_OFFSET_NU.get(int(z), RPA_Q0_OFFSET_DEFAULT)
        rel, q2vals, q2edges = self._load(fname)
        nx = rel.shape[0] - 2  # 3000

        q3bin = (q3 * 1000.0).astype(np.int64)
        q0bin = (q0 * 1000.0).astype(np.int64)
        q0bin = np.where(q0 >= 3.0, nx - 1, q0bin)
        q3bin = np.where(q3 >= 3.0, nx - 1, q3bin)
        q0bin = np.where(q0 < 0.018, 18 + offset, q0bin)
        jrow = np.clip(q0bin - offset, 0, nx + 1)
        icol = np.clip(q3bin, 0, nx + 1)
        w = rel[icol, jrow]
        w = np.where(w <= 0.001, 1.0, w)
        redo = (q0 < 0.15) & (w > 0.9)
        icol2 = np.clip(q3bin + 150, 0, nx + 1)
        w = np.where(redo, rel[icol2, jrow], w)

        q2 = q3 * q3 - q0 * q0
        w = np.where(q2 >= 9.0, 1.0, w)
        mid = (q2 > 3.0) & (q2 < 9.0)
        if mid.any():
            q2idx = np.clip(np.digitize(q2, q2edges), 0, len(q2vals) - 1)
            w = np.where(mid, q2vals[q2idx], w)
        return np.where((w >= 0.001) & (w <= 2.0), w, 1.0)

    def weight(self, q0_gev, q3_gev, mc_intType, mc_targetZ):
        q0 = np.asarray(q0_gev, dtype=np.float64)
        q3 = np.asarray(q3_gev, dtype=np.float64)
        z = np.asarray(mc_targetZ)
        out = np.ones_like(q0)
        gate = (np.asarray(mc_intType) == INT_TYPE_QE) & (z >= 6)
        for zval in np.unique(z[gate]):
            m = gate & (z == zval)
            out[m] = self._weight_one_z(q0[m], q3[m], zval)
        return out


# ---------------------------------------------------------------------------
# #4 MINOS matching efficiency — MINOSEfficiencyReweighter (ME, neutrino)
# ---------------------------------------------------------------------------
# Fifth-order polynomials in p_mu [GeV] over [1, 4] (clamped), measured at a
# low and a high beam intensity; the per-event correction is the parabola
# 1 + A*pot + B*pot^2 through those two points evaluated at the event's batch
# POT (MinosMuonEfficiencyCorrection.cxx, neutrino branch).
MINOS_POLY_LO = (0.705549, 0.499846, -0.364875, 0.131063, -0.0228281, 0.00154069)
MINOS_POLY_HI = (0.959697, -0.139194, 0.175986, -0.0802815, 0.0162377, -0.00122259)
MINOS_POT_LO = 3.9385
MINOS_POT_HI = 8.0311
MINOS_PMIN_GEV, MINOS_PMAX_GEV = 1.0, 4.0


def batch_pot(numi_pot, batch_structure, reco_vertex_batch):
    """Per-event batch POT (MinervaUniverse::GetBatchPOT, ME branch).

    k = 6 (structure 0, 3 or -1); structure 1: 4 if vertex batch < 3 else 8;
    structure 2: 5 if vertex batch < 5 else 10.
    """
    pot = np.asarray(numi_pot, dtype=np.float64)
    bs = np.asarray(batch_structure)
    vb = np.asarray(reco_vertex_batch)
    k = np.full(pot.shape, 6.0)
    k = np.where(bs == 1, np.where(vb < 3, 4.0, 8.0), k)
    k = np.where(bs == 2, np.where(vb < 5, 5.0, 10.0), k)
    return pot / k


def minos_efficiency_weight(minos_trk_p_mev, numi_pot, batch_structure,
                            reco_vertex_batch):
    """Reco-only MINOS efficiency correction (neutrino/FHC curves)."""
    p = np.clip(np.asarray(minos_trk_p_mev, dtype=np.float64) / 1000.0,
                MINOS_PMIN_GEV, MINOS_PMAX_GEV)
    corr_lo = np.polynomial.polynomial.polyval(p, MINOS_POLY_LO) - 1.0
    corr_hi = np.polynomial.polynomial.polyval(p, MINOS_POLY_HI) - 1.0
    x1, x2 = MINOS_POT_HI, MINOS_POT_LO
    det = x1 * x2 * x2 - x1 * x1 * x2
    a = (corr_hi * x2 * x2 - corr_lo * x1 * x1) / det
    b = (x1 * corr_lo - x2 * corr_hi) / det
    pot = batch_pot(numi_pot, batch_structure, reco_vertex_batch)
    return 1.0 + a * pot + b * pot * pot


# ---------------------------------------------------------------------------
# composition (Model::GetWeight at CV)
# ---------------------------------------------------------------------------
def truth_q0q3_gev(mc_incomingE, mc_primFSLepton, mc_Q2):
    """True (q0, q3) in GeV from tuple branches (TruthFunctions.h:37-47,
    BaseUniverse.cxx:77-83): q0 = Enu - Elep; q3 = sqrt(Q2 + q0^2)."""
    enu = np.asarray(mc_incomingE, dtype=np.float64)
    elep = np.asarray(mc_primFSLepton, dtype=np.float64)[:, 3]
    q0 = (enu - elep) / 1000.0
    q2 = np.asarray(mc_Q2, dtype=np.float64) / 1e6
    return q0, np.sqrt(np.maximum(q2 + q0 * q0, 0.0))


TRUTH_WEIGHT_BRANCHES = ("mc_incomingE", "mc_primFSLepton", "mc_Q2",
                         "mc_intType", "mc_targetZ",
                         "truth_genie_wgt_Rvn1pi", "truth_genie_wgt_Rvp1pi")
RECO_WEIGHT_BRANCHES = TRUTH_WEIGHT_BRANCHES + (
    "MasterAnaDev_minos_trk_p", "numi_pot", "batch_structure",
    "reco_vertex_batch")


def truth_cv_weight(arrays, flux_cv, w2p2h, rpa):
    """Product of the four truth-only CV weights (all MC fills)."""
    q0, q3 = truth_q0q3_gev(arrays["mc_incomingE"], arrays["mc_primFSLepton"],
                            arrays["mc_Q2"])
    w = flux_cv.cv_weight(np.asarray(arrays["mc_incomingE"]) / 1000.0)
    w = w * nonres_pi_weight(arrays["truth_genie_wgt_Rvn1pi"],
                             arrays["truth_genie_wgt_Rvp1pi"])
    w = w * w2p2h.weight(q0, q3, arrays["mc_intType"], arrays["mc_targetZ"])
    w = w * rpa.weight(q0, q3, arrays["mc_intType"], arrays["mc_targetZ"])
    return w


def reco_cv_weight(arrays, flux_cv, w2p2h, rpa):
    """Full CV weight for reco-side MC fills (truth product x MINOS eff)."""
    return truth_cv_weight(arrays, flux_cv, w2p2h, rpa) * minos_efficiency_weight(
        arrays["MasterAnaDev_minos_trk_p"], arrays["numi_pot"],
        arrays["batch_structure"], arrays["reco_vertex_batch"])
