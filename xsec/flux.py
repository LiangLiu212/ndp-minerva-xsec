"""Flux CV weight (MnvTune v1 weighter #1) — constrained flux / generated flux.

Implements, in Python, exactly what FluxReweighter + MnvHistoConstrainer do at
CV level (docs/cv_reweight.md §1):

1. load the PPFX-reweighted flux (``flux_E_cvweighted``) and the generated
   flux (``flux_E_unweighted``) for the playlist's flux group;
2. apply the nu-e constraint at load time (DocDB 10076 method, read from
   MnvFluxConstraint.h:180-310): per bin,
   corr = [sum_u w_u U_u / sum_u w_u] / [mean_u U_u] over the 1000 "Flux"
   error-band universes of the *_rearrangedUniverses file, with w_u from
   MParamFiles/data/FluxConstraints/sorted_NuEConstraint_FHC_RHC_IMD.txt;
   constrained CV = file CV x corr;
3. cv_weight(Enu) = constrained.Interpolate(Enu) / generated.Interpolate(Enu)
   with plain bin content above 75 GeV and weight 1 where either reads zero
   (FluxReweighter::GetFluxCVWeight);
4. integral(0, 100) with the ROOT width convention -> the normalization
   Phi_int (file units nu/m^2/POT; x1e-4 for nu/cm^2/POT).

Reading PlotUtils::MnvH1D without the MAT libraries: ROOT's
TFile::MakeProject generates + compiles the classes from the file's streamer
info once (cached next to the flux files); after gSystem.Load the objects
read normally and the universe histograms are reachable via the
fVertErrorBandMap data member. (uproot 5.7 cannot read MnvH1D — it chokes on
the empty MnvLatErrorBand map member.)

Validated: constrained Phi_int(me1D group) = 6.2299e-8 nu/cm^2/POT vs the
paper's 6.32e-8 (ratio 0.986); unconstrained would be 6.9167e-8 (+9.4%).
"""
from functools import lru_cache
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
FLUX_ROOT = REPO_ROOT / "data" / "flux"
FLUX_DIR = FLUX_ROOT / "MATFluxAndReweightFiles" / "flux"
CONSTRAINT_FILE = (FLUX_ROOT / "MParamFiles" / "data" / "FluxConstraints"
                   / "sorted_NuEConstraint_FHC_RHC_IMD.txt")
PROJ_DIR = FLUX_ROOT / "mnvh1d_proj"   # MakeProject cache (gitignored)

INTERPOLATE_MAX_GEV = 75.0  # above this, bin content (FluxReweighter.cxx slope guard)

# Playlist -> flux-file group (FluxReweighter::playlistString, ~line 1953).
FLUX_GROUP = {
    "minervame1A": "minervame1D", "minervame1B": "minervame1D",
    "minervame1C": "minervame1D", "minervame1D": "minervame1D",
    "minervame1E": "minervame1D", "minervame1F": "minervame1D",
    "minervame1G": "minervame1M", "minervame1L": "minervame1M",
    "minervame1M": "minervame1M",
    "minervame1N": "minervame1N", "minervame1O": "minervame1N",
    "minervame1P": "minervame1N",
}


def _import_root():
    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    return ROOT


def _ensure_mnvh1d_lib(ROOT, sample_file):
    """Generate (once) and load the MnvH1D classes from streamer info."""
    so = PROJ_DIR / f"{PROJ_DIR.name}.so"
    if not so.exists():
        f = ROOT.TFile.Open(str(sample_file))
        f.MakeProject(str(PROJ_DIR), "*", "recreate++")
        f.Close()
    if ROOT.gSystem.Load(str(so)) < 0:
        raise RuntimeError(f"could not load {so}")


def _read_mnvh1d(ROOT, path, hist_name, want_band=None):
    """Return (edges, cv_contents incl. under/overflow, universes or None)."""
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise FileNotFoundError(path)
    h = f.Get(hist_name)
    nb = h.GetNbinsX()
    cv = np.array([h.GetBinContent(i) for i in range(nb + 2)])
    edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)])
    unis = None
    if want_band is not None:
        bands = {str(k): v for k, v in h.fVertErrorBandMap}
        band = bands[want_band]
        unis = np.array([[band.fHists.at(u).GetBinContent(i)
                          for i in range(nb + 2)]
                         for u in range(band.fHists.size())])
    f.Close()
    return edges, cv, unis


def _load_constraint_weights(group="Flux"):
    wgts = {}
    for line in CONSTRAINT_FILE.read_text().splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("#") or parts[0] != group:
            continue
        wgts[int(parts[1])] = float(parts[2])
    return np.array([wgts[i] for i in range(len(wgts))])


class FluxCV:
    """Constrained CV flux for one playlist: weights + normalization integral."""

    def __init__(self, playlist="minervame1A"):
        group = FLUX_GROUP[playlist]
        self.playlist, self.flux_group = playlist, group
        ROOT = _import_root()
        rew_path = FLUX_DIR / f"flux-gen2thin-pdg14-{group}_rearrangedUniverses.root"
        gen_path = FLUX_DIR / f"flux-g4numiv6-pdg14-{group}.root"
        _ensure_mnvh1d_lib(ROOT, rew_path)

        self.edges, cv_rew, unis = _read_mnvh1d(ROOT, rew_path,
                                                "flux_E_cvweighted", "Flux")
        _, self.cv_generated, _ = _read_mnvh1d(ROOT, gen_path, "flux_E_unweighted")

        w = _load_constraint_weights("Flux")
        if len(w) != unis.shape[0]:
            raise ValueError(f"{len(w)} constraint weights for {unis.shape[0]} universes")
        wmean = (w[:, None] * unis).sum(axis=0) / w.sum()
        umean = unis.mean(axis=0)
        self.correction = np.divide(wmean, umean, out=np.ones_like(umean),
                                    where=umean > 0)
        self.cv_constrained = cv_rew * self.correction
        self.n_universes = unis.shape[0]

        self.cv_reweighted_unconstrained = cv_rew
        self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        self._con_real = self.cv_constrained[1:-1]
        self._gen_real = self.cv_generated[1:-1]
        self._rew_real = cv_rew[1:-1]
        # retain the raw PPFX universes + nu-e constraint weights for the flux
        # systematic (M4): the constraint enters CalcCovMx as universe weights
        # (CorrectFluxUniv -> SetUnivWgt), so the constrained flux covariance is
        # the WEIGHTED sample covariance of the raw universes.
        self._universes = unis            # (n_univ, n_bins+2), raw PPFX
        self._constraint_weights = w      # (n_univ,)
        self.constraint_weights = w       # public: SetUnivWgt for CalcCovMx

    def universe_integrals(self, emin_gev=0.0, emax_gev=100.0):
        """Per-universe flux integrals Φ_u (n_univ,) in ν/m²/POT, same width
        convention as integral()."""
        widths = np.diff(self.edges)
        lo = max(np.digitize(emin_gev, self.edges) - 1, 0)
        hi = min(np.digitize(emax_gev, self.edges) - 1, len(widths) - 1)
        real = self._universes[:, 1:-1]                      # drop under/overflow
        return (real[:, lo:hi + 1] * widths[lo:hi + 1]).sum(axis=1)

    def flux_norm_uncertainty(self, emin_gev=0.0, emax_gev=100.0, constrained=True):
        """Fractional flux-normalization uncertainty δΦ/Φ from the universe
        spread of the integrated flux. constrained=True applies the ν-e
        constraint weights (CalcCovMx weighted covariance); False = raw PPFX."""
        phi = self.universe_integrals(emin_gev, emax_gev)
        w = self._constraint_weights if constrained else np.ones_like(phi)
        mean = (w * phi).sum() / w.sum()
        var = (w * (phi - mean) ** 2).sum() / w.sum()
        return float(np.sqrt(var) / mean)

    def universe_weight_ratios(self, enu_gev):
        """(n_univ, n_events) flux-universe weight ratios to the CV flux:
        U_u(Enu) / Φ_constrained(Enu). Multiply the CV event weight by a row
        to get that universe's event weight."""
        enu = np.asarray(enu_gev, dtype=np.float64)
        cv = self._evaluate(self._con_real, enu)
        ratios = np.empty((self._universes.shape[0], enu.size))
        for u in range(self._universes.shape[0]):
            num = self._evaluate(self._universes[u, 1:-1], enu)
            ratios[u] = np.divide(num, cv, out=np.ones_like(num), where=cv != 0)
        return ratios

    def universe_ratio(self, u, enu_gev):
        """Flux weight ratio U_u(Enu)/Φ_constrained(Enu) for a SINGLE universe u
        — the memory-light per-universe form of universe_weight_ratios (the flux
        systematic pass loops 100-1000 universes over the 544k-event Truth tree,
        so the full (n_univ, n_events) matrix is avoided)."""
        enu = np.asarray(enu_gev, dtype=np.float64)
        cv = self._evaluate(self._con_real, enu)
        num = self._evaluate(self._universes[u, 1:-1], enu)
        return np.divide(num, cv, out=np.ones_like(num), where=cv != 0)

    def _evaluate(self, contents_real, enu_gev):
        """TH1::Interpolate equivalent (linear between bin centers, clamped),
        switching to plain bin content above INTERPOLATE_MAX_GEV."""
        enu = np.asarray(enu_gev, dtype=np.float64)
        interp = np.interp(enu, self.centers, contents_real)
        idx = np.clip(np.digitize(enu, self.edges) - 1, 0, len(contents_real) - 1)
        return np.where(enu > INTERPOLATE_MAX_GEV, contents_real[idx], interp)

    def cv_weight(self, enu_gev):
        """Per-event flux CV weight; 1.0 where either flux reads zero."""
        num = self._evaluate(self._con_real, enu_gev)
        den = self._evaluate(self._gen_real, enu_gev)
        ok = (num != 0) & (den != 0)
        return np.where(ok, num / np.where(ok, den, 1.0), 1.0)

    def integral(self, emin_gev=0.0, emax_gev=100.0, which="constrained"):
        """ROOT Integral(FindBin(emin), FindBin(emax), "width") over real bins.

        which: "constrained" (the normalization flux), "reweighted" (PPFX
        before the nu-e constraint, diagnostic), or "generated".
        Returns nu/m^2/POT (file units); multiply by 1e-4 for nu/cm^2/POT.
        """
        contents = {"constrained": self._con_real,
                    "reweighted": self._rew_real,
                    "generated": self._gen_real}[which]
        widths = np.diff(self.edges)
        lo = max(np.digitize(emin_gev, self.edges) - 1, 0)
        hi = min(np.digitize(emax_gev, self.edges) - 1, len(widths) - 1)
        return float((contents[lo:hi + 1] * widths[lo:hi + 1]).sum())
