"""Adapters for GiBUU FinalEvents.dat, NuHepMC (pyhepmc) and NuWro treeout."""
import tempfile
from pathlib import Path

import numpy as np
from _helpers import site, skip
from ndp.adapters.gibuu_finalevents import read_finalevents, gibuu_pdg, production_to_int_type

GIBUU_SAMPLE = """# 1:Run  2:Event  3:ID 4:Charge 5:perweight 6-8:position 9-12:momentum 13:history 14:production_ID 15:enu
      1      1    902     -1  2.0E-01  0 0 0  1.50 0.10 0.20 1.45   0  1  2.0
      1      1      1      1  0.0E+00  0 0 0  0.94 0.00 0.00 0.00   0  1  2.0
      1      1      1      1  2.0E-01  1 2 3  1.20 -0.1 -0.2 0.60   0  1  2.0
      1      2    902     -1  1.0E-01  0 0 0  0.80 0.00 0.00 0.75   0  35 2.0
      1      2      1      0  0.0E+00  0 0 0  0.94 0.00 0.00 0.00   0  35 2.0
      1      2      1      1  1.0E-01  0 0 0  1.10 0.10 0.00 0.60   0  35 2.0
      1      2    101      1  1.0E-01  0 0 0  0.30 -0.1 0.00 0.20   0  35 2.0
      2      1    902     -1  3.0E-01  0 0 0  1.00 0.00 0.00 0.90   0  34 2.0
      2      1      1      1  0.0E+00  0 0 0  0.94 0.00 0.00 0.00   0  34 2.0
"""


def test_gibuu_ids_and_process_map():
    assert gibuu_pdg(1, 1) == 2212 and gibuu_pdg(1, 0) == 2112 and gibuu_pdg(101, -1) == -211
    assert gibuu_pdg(902, -1) == 13 and gibuu_pdg(-902, 1) == -13 and gibuu_pdg(912, 0) == 14 and gibuu_pdg(999, 0) == 22
    assert [production_to_int_type(p) for p in (1, 2, 17, 32, 34, 35, 36, 37)] == [1, 2, 2, 3, 3, 5, 5, 3]


def test_gibuu_finalevents_parse():
    d = Path(tempfile.mkdtemp()); f = d / "FinalEvents.dat"; f.write_text(GIBUU_SAMPLE)
    t = read_finalevents(f, target_Z=6, target_A=12)
    assert t.n == 3 and t.has_fs
    assert list(t["int_type"]) == [1, 5, 3] and list(t["current"]) == [1, 1, 1] and list(t["nu_pdg"]) == [14, 14, 14]
    assert np.allclose(t["E_nu"], 2.0) and np.allclose(t["lep_E"], [1.5, 0.8, 1.0])
    assert np.diff(t["fs_offsets"]).tolist() == [1, 2, 0]           # struck nucleon (w=0) and lepton excluded
    assert t.meta["n_runs"] == 2 and np.allclose(t["weight"], np.array([0.2, 0.1, 0.3]) / 2)
    q2 = 2 * 2.0 * (1.5 - 1.45) - (1.5 ** 2 - (0.1 ** 2 + 0.2 ** 2 + 1.45 ** 2))
    assert abs(t["Q2"][0] - q2) < 1e-9 and t["W"][0] > 0
    assert t.norm.kind == "xsec_per_nucleon" and t.norm.xsec_per_unit_weight == 1e-38


def test_nuhepmc_roundtrip():
    try:
        import pyhepmc
    except ImportError:
        skip("pyhepmc not installed in this interpreter (use pixi run test)")
    from ndp.adapters.nuhepmc import read_nuhepmc
    d = Path(tempfile.mkdtemp()); f = d / "toy.hepmc"
    ri = pyhepmc.GenRunInfo()
    ri.attributes["NuHepMC.ProcessInfo[251].Name"] = pyhepmc.StringAttribute("QESpectralCC1p0pi")
    ri.attributes["NuHepMC.ProcessInfo[452].Name"] = pyhepmc.StringAttribute("RES_Spectral_Func")
    ri.attributes["NuHepMC.Units.CrossSection.Unit"] = pyhepmc.StringAttribute("pb")
    ri.attributes["NuHepMC.Units.CrossSection.TargetScale"] = pyhepmc.StringAttribute("PerAtom")
    ri.attributes["NuHepMC.FluxAveragedTotalCrossSection"] = pyhepmc.DoubleAttribute(1.2e5)   # pb / C12 atom
    ri.weight_names = ["CV"]
    with pyhepmc.open(str(f), "w") as out:
        for pid, lep_pdg in ((251, 13), (452, 13)):
            evt = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM); evt.run_info = ri
            v = pyhepmc.GenVertex()
            nu = pyhepmc.GenParticle((0, 0, 2.0, 2.0), 14, 4); tgt = pyhepmc.GenParticle((0, 0, 0, 11.18), 1000060120, 20)
            mu = pyhepmc.GenParticle((0.3, 0.0, 1.3, 1.34), lep_pdg, 1); pr = pyhepmc.GenParticle((-0.3, 0.0, 0.7, 1.2), 2212, 1)
            for p in (nu, tgt):
                v.add_particle_in(p)
            for p in (mu, pr):
                v.add_particle_out(p)
            evt.add_vertex(v); evt.weights = [1.0]
            evt.attributes["signal_process_id"] = pyhepmc.IntAttribute(pid)
            out.write(evt)
    t = read_nuhepmc(f)
    assert t.n == 2 and list(t["int_type"]) == [1, 2] and list(t["nu_pdg"]) == [14, 14] and list(t["target_A"]) == [12, 12]
    assert np.allclose(t["lep_E"], 1.34) and np.diff(t["fs_offsets"]).tolist() == [1, 1]
    assert t.norm.kind == "xsec_per_nucleon"
    assert abs(t.norm.xsec_per_unit_weight * t["weight"].sum() - 1.2e5 * 1e-36 / 12) < 1e-45


def test_nuwro_smoke_output_if_present():
    p = site().repo_root / "external/generator-smoke/nuwro/smoke.root"
    if not p.exists():
        skip("no NuWro smoke output (pixi run test-generators)")
    try:
        import ROOT  # noqa: F401
    except ImportError:
        skip("PyROOT not in this interpreter (use pixi run test)")
    import os
    os.environ.setdefault("NUWRO", str(site().repo_root / "external/nuwro"))
    from ndp.adapters.nuwro_root import read_nuwro
    t = read_nuwro(p, entry_stop=50)
    assert t.n == 50 and set(t["nu_pdg"]) == {14} and t.norm.kind == "xsec_per_nucleon" and t.has_fs
