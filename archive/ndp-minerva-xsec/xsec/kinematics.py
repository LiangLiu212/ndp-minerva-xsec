"""Vectorized muon kinematics and detector geometry.

NumPy arrays in, NumPy arrays out; no I/O. Formulas are lifted from the
certified/frozen exploration-repo code so results are bit-compatible with the
golden counts:

- theta3d:        tools/cc_inclusive_selector.py (certified, Stage-B parity)
- hexagon test:   tools/cc_inclusive_selector.py apothem_pass
- beam rotation + true theta/p: exploring/dsigma_dpt.py true_theta_and_p
                  (RotateX convention, TruthFunctions.h:21-27)
"""
import numpy as np

from .constants import (APOTHEM_INTERCEPT_MM, APOTHEM_MM, APOTHEM_SLOPE,
                        NUMI_BEAM_ANGLE_RAD)


def theta3d(theta_x, theta_y):
    """3D polar angle from the two projected track angles (radians).

    acos(1/sqrt(1 + tan^2(tx) + tan^2(ty))) — selector's theta3d, vectorized.

    Cross-checked vs MAT BaseUniverse::theta3D (BaseUniverse.cxx:67-75,
    called by GetThetamu via MuonFunctions.h:133): MAT uses
    acos(1/sqrt(sec^2(tx) + sec^2(ty) - 1)), algebraically IDENTICAL since
    sec^2 - 1 = tan^2. MAT additionally folds backward tracks
    (theta_x or theta_y > pi/2 -> pi - angle); the certified selector and
    this port deliberately omit the fold — a backward muon cannot satisfy
    the MINOS match, so the final selection is unaffected (certified parity:
    0 mismatches on the golden files).
    """
    tx = np.tan(np.asarray(theta_x, dtype=np.float64))
    ty = np.tan(np.asarray(theta_y, dtype=np.float64))
    return np.arccos(1.0 / np.sqrt(1.0 + tx * tx + ty * ty))


def in_hex_apothem(x, y):
    """Flat-top hexagonal fiducial test (apothem 850 mm), strict <.

    abs(x) < apothem  AND  abs(y) < slope*abs(x) + intercept
    """
    ax = np.abs(np.asarray(x, dtype=np.float64))
    ay = np.abs(np.asarray(y, dtype=np.float64))
    return (ax < APOTHEM_MM) & (ay < APOTHEM_SLOPE * ax + APOTHEM_INTERCEPT_MM)


def reco_pt_pz_gev(leptonE, theta_x, theta_y):
    """Reconstructed muon (p_T, p_parallel) in GeV/c — Stage-B convention.

    p = |leptonE[:, 0:3]| (MeV); theta = theta3d(muon_thetaX, muon_thetaY);
    p_T = p*sin(theta)/1000, p_parallel = p*cos(theta)/1000.
    Source: dsigma_dpt.py:179-190 (frozen; pT identical to Stage B
    stageB_pT_selection.py:189-194); p_parallel uses the same |p| and angle.
    """
    lep = np.asarray(leptonE, dtype=np.float64)
    p_mev = np.sqrt(lep[:, 0] ** 2 + lep[:, 1] ** 2 + lep[:, 2] ** 2)
    theta = theta3d(theta_x, theta_y)
    return (p_mev * np.sin(theta) / 1000.0,
            p_mev * np.cos(theta) / 1000.0)


def true_theta_p(px, py, pz):
    """True angle w.r.t. the NuMI beam axis and momentum magnitude (MeV/c).

    Beam-frame z component: zb = sin(a)*py + cos(a)*pz with
    a = NUMI_BEAM_ANGLE_RAD (RotateX convention); theta = acos(zb/p), clipped.
    Rows with p == 0 return (theta=0, p=0), matching the scalar frozen code
    (they subsequently fail the p_z phase-space cut).
    """
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64)
    pz = np.asarray(pz, dtype=np.float64)
    c, s = np.cos(NUMI_BEAM_ANGLE_RAD), np.sin(NUMI_BEAM_ANGLE_RAD)
    zb = s * py + c * pz
    p = np.sqrt(px * px + py * py + pz * pz)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(p > 0.0, zb / np.where(p > 0.0, p, 1.0), 1.0)
    theta = np.arccos(np.clip(ratio, -1.0, 1.0))
    theta = np.where(p > 0.0, theta, 0.0)
    return theta, np.where(p > 0.0, p, 0.0)
