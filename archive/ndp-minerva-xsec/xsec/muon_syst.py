"""Per-event muon-reconstruction systematic shifts (lateral universes).

Ports of the MAT-MINERvA muon-systematics universe classes. The muon
reconstruction category (Fig 8 "Muon Reconstruction") is the dominant
systematic; it is the quadrature sum of several ±1σ bands, each a reco-side
shift re-run through the full chain. Magnitudes are from NSFDefaults.h.

Two kinds of shift:
  * momentum scale — |p| changes, the track angle is fixed, so reco p_T and
    p_∥ both scale by the same factor p_new/p_total (selection is angle-only,
    so it is unchanged → re-bin only);
  * track angle — the angle changes at fixed |p|, so reco p_T/p_∥ are
    recomputed AND the 20° muon-angle cut must be re-applied (re-select).

The total reco momentum splits p_total = p_minerva + p_minos, with
p_minos = MasterAnaDev_minos_trk_p and p_total = |MasterAnaDev_leptonE|
(MuonFunctions.h:57-112).
"""
import numpy as np

from .constants import NUMI_BEAM_ANGLE_RAD

# --- NSFDefaults.h magnitudes ----------------------------------------------
# Muon_Energy_MINERvA: absolute MeV shift on the MINERvA momentum component.
# Our fiducial (vtx Z 5980-8422) is entirely past nuclearTargetZEnd=5835, so
# every event uses the "NoNuke" dE/dx (30) + material-assay (11) values
# (MuonSystematics.cxx:81-91).
DEDX_UNCERT_NONUKE_MEV = 30.0
MATERIAL_ASSAY_NONUKE_MEV = 11.0
MUON_ENERGY_MINERVA_SHIFT_MEV = float(
    np.hypot(DEDX_UNCERT_NONUKE_MEV, MATERIAL_ASSAY_NONUKE_MEV))  # 31.95 MeV

# Muon_Energy_MINOS: fractional shift on the MINOS component (range always,
# curvature only if reconstructed by curvature; high/low-p split at 1 GeV).
MINOS_RANGE_ERR = 0.00984
MINOS_HIGHP_CURV_ERR = 0.006
MINOS_LOWP_CURV_ERR = 0.025
MINOS_CURV_PTHRESH_MEV = 1000.0

MUON_RESOLUTION_ERR = 0.004          # fractional, on (p_true - p_reco)
BEAM_THETAX_ERR = 0.001              # rad
BEAM_THETAY_ERR = 0.0009             # rad
MUON_ANGLE_FRAC_RES = 0.02           # fractional, on (theta_true - theta_reco)


# --- momentum-scale bands (return p_new / p_total) -------------------------
def pmu_minerva_scale(p_total_mev, nsigma):
    """Muon_Energy_MINERvA: total p shifts by ±31.95 MeV absolute
    (MuonSystematics.cxx:89-103). Momentum-dependent fractional scale."""
    p = np.asarray(p_total_mev, np.float64)
    return 1.0 + nsigma * MUON_ENERGY_MINERVA_SHIFT_MEV / p


def pmu_minos_scale(p_total_mev, p_minos_mev, used_curvature, nsigma):
    """Muon_Energy_MINOS: the MINOS component shifts fractionally by
    f=√(range²+curv²) (MuonSystematics.cxx:122-152); curv applies only to
    curvature-reconstructed tracks. Total scale = 1 + nσ·f·p_minos/p_total."""
    p = np.asarray(p_total_mev, np.float64)
    pm = np.asarray(p_minos_mev, np.float64)
    curv = np.where(np.asarray(used_curvature).astype(bool),
                    np.where(pm > MINOS_CURV_PTHRESH_MEV,
                             MINOS_HIGHP_CURV_ERR, MINOS_LOWP_CURV_ERR), 0.0)
    f = np.sqrt(MINOS_RANGE_ERR ** 2 + curv ** 2)
    return 1.0 + nsigma * f * pm / p


def pmu_resolution_scale(p_total_mev, p_true_mev, nsigma):
    """Muon momentum resolution (MuonResolutionSystematics.cxx:65-87):
    p shifts by (p_true − p_reco)·nσ·0.004."""
    p = np.asarray(p_total_mev, np.float64)
    pt = np.asarray(p_true_mev, np.float64)
    return 1.0 + (pt - p) * (nsigma * MUON_RESOLUTION_ERR) / p


# --- track-angle bands (return shifted theta_x, theta_y) -------------------
def beam_angle_shift(theta_x, theta_y, axis, nsigma):
    """BeamAngleX/Y (AngleSystematics.cxx:48-74): add nσ·err to one projection
    (X err 0.001, Y err 0.0009 rad)."""
    tx = np.array(theta_x, np.float64, copy=True)
    ty = np.array(theta_y, np.float64, copy=True)
    if axis == "x":
        tx += nsigma * BEAM_THETAX_ERR
    else:
        ty += nsigma * BEAM_THETAY_ERR
    return tx, ty


def true_theta_xy(mc_primFSLepton):
    """True track angles (θ_x, θ_y) in the beam frame
    (MuonResolutionSystematics.cxx GetTrueTheta{X,Y}mu): RotateX by the NuMI
    beam angle, then per-projection signed acos(z'/√(proj²+z'²))."""
    lep = np.asarray(mc_primFSLepton, np.float64)
    px, py, pz = lep[:, 0], lep[:, 1], lep[:, 2]
    c, s = np.cos(NUMI_BEAM_ANGLE_RAD), np.sin(NUMI_BEAM_ANGLE_RAD)
    xr = px                       # RotateX leaves x unchanged
    yr = py * c - pz * s
    zr = py * s + pz * c          # == frozen true_theta_p z-component
    dx = np.sqrt(xr * xr + zr * zr)
    dy = np.sqrt(yr * yr + zr * zr)
    tx = np.sign(xr) * np.arccos(np.clip(zr / np.where(dx > 0, dx, 1.0), -1, 1))
    ty = np.sign(yr) * np.arccos(np.clip(zr / np.where(dy > 0, dy, 1.0), -1, 1))
    return tx, ty


def angle_resolution_shift(theta_x, theta_y, true_x, true_y, axis, nsigma):
    """MuonAngle{X,Y}Resolution (MuonResolutionSystematics.cxx:113-169):
    θ shifts by (θ_true − θ_reco)·nσ·0.02 on the chosen projection."""
    tx = np.array(theta_x, np.float64, copy=True)
    ty = np.array(theta_y, np.float64, copy=True)
    if axis == "x":
        tx += (np.asarray(true_x, np.float64) - tx) * (nsigma * MUON_ANGLE_FRAC_RES)
    else:
        ty += (np.asarray(true_y, np.float64) - ty) * (nsigma * MUON_ANGLE_FRAC_RES)
    return tx, ty
