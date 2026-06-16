"""Tests for xsec.weights (CV weights #2-#5).

Synthetic hand-traces against the C++ formulas; RPA/2p2h need the unpacked
tarball (skipped otherwise). The batch-POT units check streams the golden
data file (network marker).
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from xsec import weights
from xsec.constants import NONRES_PI_WEIGHT

REPO_ROOT = Path(__file__).resolve().parents[1]

needs_files = pytest.mark.skipif(not weights.REWEIGHT_DIR.exists(),
                                 reason="tarball not unpacked under data/flux/")


# ------------------------------------------------------------ non-res pi (#2)
def test_nonres_pi_truth_table():
    seven = lambda v2: [1.0, 1.0, v2, 1.0, 1.0, 1.0, 1.0]
    rvn = np.array([seven(0.6), seven(1.0), seven(1.0), seven(0.99)])
    rvp = np.array([seven(1.0), seven(0.4), seven(1.0), seven(1.0)])
    w = weights.nonres_pi_weight(rvn, rvp)
    assert w.tolist() == [NONRES_PI_WEIGHT, NONRES_PI_WEIGHT, 1.0, NONRES_PI_WEIGHT]


def test_nonres_pi_variation_ratio():
    seven = lambda v2: [1.0, 1.0, v2, 1.0, 1.0, 1.0, 1.0]
    rvn = np.array([seven(0.6), seven(1.0), seven(1.0)])   # tagged, untagged, untagged
    rvp = np.array([seven(1.0), seven(0.4), seven(1.0)])   # untagged, tagged, untagged
    f = 0.04 / 0.43
    rp = weights.nonres_pi_variation_ratio(rvn, rvp, +1.)
    rm = weights.nonres_pi_variation_ratio(rvn, rvp, -1.)
    assert np.allclose(rp, [1 + f, 1 + f, 1.0])
    assert np.allclose(rm, [1 - f, 1 - f, 1.0])
    # CV × ratio reproduces the 0.43 ± 0.04 weight on tagged events
    assert np.isclose(NONRES_PI_WEIGHT * rp[0], 0.47)
    assert np.isclose(NONRES_PI_WEIGHT * rm[0], 0.39)


# ----------------------------------------------------------------- 2p2h (#3)
@needs_files
def test_2p2h_params_parsed():
    w = weights.TwoP2HWeight()
    assert (round(w.norm, 4), round(w.meanq0, 6), round(w.meanq3, 5)) == \
           (10.5798, 0.254032, 0.50834)
    assert (round(w.sigmaq0, 7), round(w.sigmaq3, 6), round(w.corr, 6)) == \
           (0.0571035, 0.129051, 0.875287)


@needs_files
def test_2p2h_weight_values_and_gating():
    w = weights.TwoP2HWeight()
    # at the Gaussian mean: w = 1 + norm, MEC on carbon only
    q0 = np.array([w.meanq0] * 4)
    q3 = np.array([w.meanq3] * 4)
    int_t = np.array([8, 1, 8, 8])
    z = np.array([6, 6, 1, 6])
    out = w.weight(q0, q3, int_t, z)
    assert abs(out[0] - (1.0 + w.norm)) < 1e-12
    assert out[1] == 1.0          # QE not weighted at CV
    assert out[2] == 1.0          # hydrogen excluded
    assert abs(out[3] - out[0]) < 1e-12
    # far tail -> ~1
    far = w.weight(np.array([2.5]), np.array([2.9]), np.array([8]), np.array([6]))
    assert abs(far[0] - 1.0) < 1e-6
    # hand-traced Gaussian at an off-center point
    q0p, q3p = 0.30, 0.60
    rho, s0, s3 = w.corr, w.sigmaq0, w.sigmaq3
    zarg = ((q0p - w.meanq0) ** 2 / s0 ** 2 + (q3p - w.meanq3) ** 2 / s3 ** 2
            - 2 * rho * (q0p - w.meanq0) * (q3p - w.meanq3) / (s0 * s3))
    expect = 1.0 + w.norm * math.exp(-0.5 * zarg / (1 - rho ** 2))
    got = w.weight(np.array([q0p]), np.array([q3p]), np.array([8]), np.array([6]))[0]
    assert abs(got - expect) < 1e-12


# ------------------------------------------------------------------ RPA (#5)
@needs_files
def test_rpa_gating_and_traps():
    rpa = weights.RPAWeight()
    # non-QE and light targets are exactly 1
    out = rpa.weight(np.array([0.1, 0.1]), np.array([0.4, 0.4]),
                     np.array([2, 1]), np.array([6, 1]))
    assert out.tolist() == [1.0, 1.0]
    # Q2 >= 9 -> 1 even for QE on carbon
    out = rpa.weight(np.array([0.5]), np.array([3.05]), np.array([1]), np.array([6]))
    assert out[0] == 1.0
    # typical low-Q2 QE on carbon: suppression (w < 1), within sanity range
    out = rpa.weight(np.array([0.08]), np.array([0.35]), np.array([1]), np.array([6]))
    assert 0.001 <= out[0] <= 2.0 and out[0] < 1.0


@needs_files
def test_rpa_matches_manual_lookup():
    """Replicate getWeightInternal by hand on points covering each branch."""
    import uproot
    rpa = weights.RPAWeight()
    f = uproot.open(weights.REWEIGHT_DIR / "outNievesRPAratio-nu12C-20GeV-20170202.root")
    rel = f["hrelratio"].values(flow=True)
    q2v = f["hQ2relratio"].values(flow=True)
    q2e = f["hQ2relratio"].axis(0).edges()
    off = weights.RPA_Q0_OFFSET_NU[6]

    def manual(q0, q3):
        nx = 3000
        q3b = int(q3 * 1000) if q3 < 3.0 else nx - 1
        q0b = int(q0 * 1000) if q0 < 3.0 else nx - 1
        if q0 < 0.018:
            q0b = 18 + off
        w = rel[min(q3b, nx + 1), min(max(q0b - off, 0), nx + 1)]
        if w <= 0.001:
            w = 1.0
        if q0 < 0.15 and w > 0.9:
            w = rel[min(q3b + 150, nx + 1), min(max(q0b - off, 0), nx + 1)]
        q2 = q3 * q3 - q0 * q0
        if q2 >= 9.0:
            w = 1.0
        elif q2 > 3.0:
            w = q2v[min(np.digitize(q2, q2e), len(q2v) - 1)]
        return w if 0.001 <= w <= 2.0 else 1.0

    pts = [(0.05, 0.3), (0.01, 0.25), (0.12, 0.5), (0.4, 0.9), (1.2, 2.2),
           (0.5, 2.0), (0.3, 2.5), (2.0, 2.9), (3.5, 3.5), (0.08, 0.35)]
    q0s = np.array([p[0] for p in pts])
    q3s = np.array([p[1] for p in pts])
    got = rpa.weight(q0s, q3s, np.full(len(pts), 1), np.full(len(pts), 6))
    for (q0, q3), g in zip(pts, got):
        assert abs(g - manual(q0, q3)) < 1e-12, (q0, q3, g, manual(q0, q3))


@needs_files
def test_rpa_variation_bands_manual_lookup():
    """Parity for the 4 RPA error bands vs a hand replication of
    weightRPA.cxx:134-235 (suppression ±25 %, enhancement via the non-rel ratio)."""
    import uproot
    rpa = weights.RPAWeight()
    f = uproot.open(weights.REWEIGHT_DIR / "outNievesRPAratio-nu12C-20GeV-20170202.root")
    rel = f["hrelratio"].values(flow=True)
    nonrel = f["hnonrelratio"].values(flow=True)
    q2rel = f["hQ2relratio"].values(flow=True)
    q2nonrel = f["hQ2nonrelratio"].values(flow=True)
    q2e = f["hQ2relratio"].axis(0).edges()
    off = weights.RPA_Q0_OFFSET_NU[6]
    nx = 3000

    def look(hist, q1, q0, q3):
        q3b = int(q3 * 1000) if q3 < 3.0 else nx - 1
        q0b = (18 + off) if q0 < 0.018 else (int(q0 * 1000) if q0 < 3.0 else nx - 1)
        jrow = min(max(q0b - off, 0), nx + 1)
        w = hist[min(q3b, nx + 1), jrow]
        if w <= 0.001:
            w = 1.0
        if q0 < 0.15 and w > 0.9:
            w = hist[min(q3b + 150, nx + 1), jrow]
        q2 = q3 * q3 - q0 * q0
        if q2 >= 9.0:
            w = 1.0
        elif q2 > 3.0:
            w = q1[min(np.digitize(q2, q2e), len(q1) - 1)]
        return w if 0.001 <= w <= 2.0 else 1.0

    pts = [(0.05, 0.3), (0.12, 0.5), (0.08, 0.35), (0.30, 0.90), (0.50, 1.20)]
    for q0, q3 in pts:
        cv = look(rel, q2rel, q0, q3)
        ext = look(nonrel, q2nonrel, q0, q3)
        q2 = q3 * q3 - q0 * q0
        elp = cv + 0.25 * (1 - cv) if cv < 1.0 else cv
        elm = cv - 0.25 * (1 - cv) if cv < 1.0 else cv
        ep = cv + 0.6 * (ext - cv)
        if q2 < 0.9:
            ep += 1.5 * (0.9 - q2) * (ext - ep)
        ep = min(ep, ext); ep = max(ep, cv + 0.03)
        em = cv - 0.6 * (ext - cv); em = min(em, cv - 0.03)
        if q2 > 1.0 and em < 1.0:
            em = 1.0
        g = rpa._weights_one_z(np.array([q0]), np.array([q3]), 6)
        for k, exp in enumerate((cv, elp, elm, ep, em)):
            assert abs(g[k][0] - exp) < 1e-9, (q0, q3, k, g[k][0], exp)

    # CV from the band routine equals the standalone CV path
    q0s = np.array([p[0] for p in pts]); q3s = np.array([p[1] for p in pts])
    cvs = rpa._weights_one_z(q0s, q3s, 6)[0]
    assert np.allclose(cvs, rpa._weight_one_z(q0s, q3s, 6), atol=1e-12)
    # ratio is 1 outside the QE/Z>=6 gate
    r = rpa.variation_ratio(np.array([0.1]), np.array([0.4]),
                            np.array([2]), np.array([6]), "HighQ2", 1)
    assert r[0] == 1.0


# ------------------------------------------------------------ MINOS eff (#4)
@needs_files
def test_twop2h_variation_gating_and_files():
    q0 = np.array([0.254, 0.254, 0.254, 0.254])
    q3 = np.array([0.508, 0.508, 0.508, 0.508])
    it = np.array([8, 8, 8, 1])                 # MEC, MEC, MEC, QE
    z = np.array([6, 6, 6, 6])
    nuc = np.array([2000000200, 2000000201, 2000000202, 2000000200])  # nn, np, pp, -
    w1 = weights.twop2h_variation_weight(q0, q3, it, z, nuc, 1)   # nn/pp only
    assert w1[0] != 1.0 and w1[2] != 1.0 and w1[1] == 1.0 and w1[3] == 1.0
    w2 = weights.twop2h_variation_weight(q0, q3, it, z, nuc, 2)   # np only
    assert w2[1] != 1.0 and w2[0] == 1.0
    w3 = weights.twop2h_variation_weight(q0, q3, it, z, nuc, 3)   # QE->2p2h
    assert w3[3] != 1.0 and w3[0] == 1.0        # acts on QE row only


@needs_files
def test_twop2h_variation_ratio_swaps_cv():
    """The universe weight is CV_stack × (variation / CV) — the ratio swaps the
    CV 2p2h factor for the mode's variation (MnvTuneSystematics.cxx:18-61)."""
    w2p2h = weights.TwoP2HWeight()
    q0 = np.array([0.254, 0.254, 0.254])
    q3 = np.array([0.508, 0.508, 0.508])
    it = np.array([8, 8, 1])                              # MEC nn, MEC np, QE
    z = np.array([6, 6, 6])
    nuc = np.array([2000000200, 2000000201, 2000000200])  # nn, np, -
    cv = w2p2h.weight(q0, q3, it, z)
    for mode in (1, 2, 3):
        var = weights.twop2h_variation_weight(q0, q3, it, z, nuc, mode)
        ratio = weights.twop2h_variation_ratio(q0, q3, it, z, nuc, w2p2h, mode)
        assert np.allclose(ratio, var / cv)
    # QE row: CV 2p2h is 1, so the mode-3 ratio is the bare QE variation (>1);
    # an MEC row under mode 3 has variation off (1.0) -> ratio 1/cv < 1.
    r3 = weights.twop2h_variation_ratio(q0, q3, it, z, nuc, w2p2h, 3)
    assert cv[2] == 1.0 and r3[2] > 1.0
    assert r3[0] < 1.0


def test_minos_efficiency_error_and_universe():
    # theta-dependent fractional error, ~1-3%, rises with angle
    err0 = weights.minos_efficiency_error(np.array([0.0]))[0]
    err20 = weights.minos_efficiency_error(np.array([20.0]))[0]
    assert 0.005 < err0 < 0.05 and err20 > err0
    # clamp beyond 40 deg
    assert weights.minos_efficiency_error(np.array([50.0]))[0] == \
           weights.minos_efficiency_error(np.array([40.0]))[0]
    # ±sigma universe brackets the CV
    args = (np.array([3000.0]), np.array([180.0]), np.array([0]), np.array([0]))
    cv = weights.minos_efficiency_weight(*args)
    up = weights.minos_efficiency_universe(*args, np.array([10.0]), +1)
    dn = weights.minos_efficiency_universe(*args, np.array([10.0]), -1)
    assert dn[0] < cv[0] < up[0]
    assert np.isclose((up[0] + dn[0]) / 2, cv[0])


def test_batch_pot_divisors():
    pot = np.full(6, 60.0)
    bs = np.array([0, 3, -1, 1, 1, 2])
    vb = np.array([0, 0, 0, 2, 3, 5])
    expect = [10.0, 10.0, 10.0, 15.0, 7.5, 6.0]
    assert weights.batch_pot(pot, bs, vb).tolist() == expect


def test_minos_correction_reproduces_anchor_curves():
    """By construction the parabola passes through (pot_lo, corr_lo(p)) and
    (pot_hi, corr_hi(p)) — check both anchors at several momenta."""
    for p_gev in (1.0, 1.7, 2.5, 3.3, 4.0):
        lo = np.polynomial.polynomial.polyval(p_gev, weights.MINOS_POLY_LO)
        hi = np.polynomial.polynomial.polyval(p_gev, weights.MINOS_POLY_HI)
        for pot, expect in ((weights.MINOS_POT_LO, lo), (weights.MINOS_POT_HI, hi)):
            got = weights.minos_efficiency_weight(
                np.array([p_gev * 1000.0]), np.array([pot * 6.0]),
                np.array([0]), np.array([0]))[0]
            assert abs(got - expect) < 1e-10, (p_gev, pot, got, expect)


def test_minos_momentum_clamp():
    w_low = weights.minos_efficiency_weight(np.array([500.0]), np.array([30.0]),
                                            np.array([0]), np.array([0]))
    w_min = weights.minos_efficiency_weight(np.array([1000.0]), np.array([30.0]),
                                            np.array([0]), np.array([0]))
    assert w_low[0] == w_min[0]


def test_genie_knob_ratio():
    # 7-elem knob [-3,-2,-1,CV,+1,+2,+3] sigma; CV at idx 3
    knob = np.array([[0.90, 0.92, 0.95, 1.0, 1.06, 1.12, 1.20],
                     [1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0]])
    assert np.allclose(weights.genie_knob_ratio(knob, +1), [1.06, 1.5])
    assert np.allclose(weights.genie_knob_ratio(knob, -1), [0.95, 0.5])
    # CV entry 0 -> ratio 1 (guarded)
    zero = np.array([[1.0, 1.0, 1.0, 0.0, 2.0, 1.0, 1.0]])
    assert weights.genie_knob_ratio(zero, +1)[0] == 1.0


def test_genie_knob_ratio_jagged():
    # uproot may hand back an object (jagged) array; _knob_element handles it
    knob = np.empty(2, dtype=object)
    knob[0] = np.array([0.9, 0.92, 0.95, 1.0, 1.06, 1.12, 1.2])
    knob[1] = np.array([1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    assert np.allclose(weights.genie_knob_ratio(knob, +1), [1.06, 1.5])


@pytest.mark.network
def test_batch_pot_units_in_real_data():
    """Stream the golden data file's beam branches: the batch POT must land
    in the intensity range the correction curves were measured at (units of
    1e12 POT: curves anchored at 3.94 and 8.03)."""
    import uproot
    spec = json.loads((REPO_ROOT / "config" / "datasets" / "me1A_single_pair.json").read_text())
    url = next(e["url"] for e in spec["files"] if e["role"] == "data")
    try:
        f = uproot.open(url)
    except Exception as err:
        pytest.skip(f"xrootd unreachable: {err}")
    with f:
        arrs = f["MasterAnaDev"].arrays(["numi_pot", "batch_structure",
                                         "reco_vertex_batch"], library="np")
    bp = weights.batch_pot(arrs["numi_pot"], arrs["batch_structure"],
                           arrs["reco_vertex_batch"])
    med = float(np.median(bp))
    assert 2.0 < med < 12.0, f"median batch POT {med} — units assumption broken"
