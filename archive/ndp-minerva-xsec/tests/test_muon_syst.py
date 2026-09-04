"""Tests for xsec.muon_syst — per-event muon-reconstruction shifts.

Hand-traces against the MAT-MINERvA universe formulas (MuonSystematics.cxx,
AngleSystematics.cxx, MuonResolutionSystematics.cxx) + NSFDefaults.h.
"""
import numpy as np

from xsec import muon_syst as ms
from xsec.constants import NUMI_BEAM_ANGLE_RAD


def test_minerva_shift_is_absolute_3195_mev():
    assert abs(ms.MUON_ENERGY_MINERVA_SHIFT_MEV - np.hypot(30.0, 11.0)) < 1e-12
    p = np.array([5000.0, 2000.0])
    up = ms.pmu_minerva_scale(p, +1.)
    dn = ms.pmu_minerva_scale(p, -1.)
    # absolute MeV shift -> momentum-dependent fractional, symmetric about 1
    assert np.allclose((up - 1.) * p, ms.MUON_ENERGY_MINERVA_SHIFT_MEV)
    assert np.allclose((up + dn) / 2, 1.0)


def test_minos_scale_range_only_vs_curvature():
    p = np.array([5000.0, 5000.0, 5000.0])
    pm = np.array([4000.0, 4000.0, 500.0])
    curv = np.array([False, True, True])
    up = ms.pmu_minos_scale(p, pm, curv, +1.)
    # range-only: f = 0.00984
    assert np.isclose(up[0], 1 + 0.00984 * 4000 / 5000)
    # high-p curvature: f = sqrt(0.00984^2 + 0.006^2)
    f_hi = np.hypot(0.00984, 0.006)
    assert np.isclose(up[1], 1 + f_hi * 4000 / 5000)
    # low-p curvature (<1 GeV): f = sqrt(0.00984^2 + 0.025^2)
    f_lo = np.hypot(0.00984, 0.025)
    assert np.isclose(up[2], 1 + f_lo * 500 / 5000)


def test_resolution_scale_toward_truth():
    p = np.array([5000.0])
    ptrue = np.array([5100.0])
    up = ms.pmu_resolution_scale(p, ptrue, +1.)
    assert np.isclose(up[0], 1 + (5100 - 5000) * 0.004 / 5000)
    # nsigma=0 is the identity
    assert np.isclose(ms.pmu_resolution_scale(p, ptrue, 0.)[0], 1.0)


def test_beam_angle_shift_one_projection():
    tx = np.array([0.05, 0.10]); ty = np.array([0.02, 0.03])
    sx, sy = ms.beam_angle_shift(tx, ty, "x", +1.)
    assert np.allclose(sx, tx + 0.001) and np.allclose(sy, ty)
    sx, sy = ms.beam_angle_shift(tx, ty, "y", -1.)
    assert np.allclose(sx, tx) and np.allclose(sy, ty - 0.0009)


def test_true_theta_xy_matches_rotatex_definition():
    # one muon, beam-frame angles via explicit RotateX
    lep = np.array([[300.0, 120.0, 5000.0, 5010.0]])
    tx, ty = ms.true_theta_xy(lep)
    c, s = np.cos(NUMI_BEAM_ANGLE_RAD), np.sin(NUMI_BEAM_ANGLE_RAD)
    xr, yr, zr = 300.0, 120.0 * c - 5000.0 * s, 120.0 * s + 5000.0 * c
    ex = np.sign(xr) * np.arccos(zr / np.hypot(xr, zr))
    ey = np.sign(yr) * np.arccos(zr / np.hypot(yr, zr))
    assert np.isclose(tx[0], ex) and np.isclose(ty[0], ey)


def test_angle_resolution_shift_toward_truth():
    tx = np.array([0.050]); ty = np.array([0.020])
    truex = np.array([0.052]); truey = np.array([0.018])
    sx, sy = ms.angle_resolution_shift(tx, ty, truex, truey, "x", +1.)
    assert np.isclose(sx[0], 0.050 + (0.052 - 0.050) * 0.02) and np.allclose(sy, ty)
