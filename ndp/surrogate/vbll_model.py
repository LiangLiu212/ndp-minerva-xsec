"""Inference-side port of the VBLL_SurrogateModel `ParticleSurrogate`: (particle type, true 4-vector)
-> reconstructed 4-vector with a per-event sigma, loadable from the platform's own `arrays.npz` +
`surrogate.json` without the VBLL repository on the path.

Source model: github.com/fnal-nd-aiml/VBLL_SurrogateModel, branch sdey2_dev @ 77c4764, `code/model.py`.
    Embedding(2, d_embed) of the particle type (0 muon, 1 proton) ++ normalised truth (E, px, py, pz)
    -> MLP(hidden x n_layers, ReLU) = h
    -> one vbll head per particle: HetRegression (head_type 'het') or Regression ('standard')
    (+ a DiscClassification head on h that predicts the particle type; unused for smearing, kept so the
       checkpoint loads strictly and the architecture is provably the trained one).
Units: the network sees the training CSV's *normalised* units; `normaliser` (train-split mean / std per
particle and component, inputs = truth, outputs = reco) converts MeV <-> normalised. The frame is whatever
the training CSV used (recorded in the spec; pinned by the plan's Step 2).

Predictive distribution, vbll 0.4.9 `HetRegression.predictive_sample(x, consistent_variance=False)`:
    mean  = W.mean @ h
    sigma = sqrt((Var[W @ h] + 1) * exp(s)),   s ~ Normal(M.mean @ h, sqrt(exp(M_logdiag) . h^2))
i.e. the noise map M is a posterior and sigma is *drawn* per call: the reference evaluation
(`evaluate.collect_predictions`) draws it once per batch, so its CV / coverage / outlier numbers are
stochastic; the mean is deterministic. The platform smearer draws a fresh sigma for every sample.
For the standard head, sigma = sqrt(Var[W @ h] + exp(noise_logdiag)) is deterministic.

torch and vbll are imported lazily so the default (torch-free) environment can still import ndp.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

PARTICLES = ("muon", "proton")
COMPONENTS = ("E", "px", "py", "pz")
PARTICLE_INDEX = {"muon": 0, "proton": 1}
KIND = "vbll_model"


@dataclass
class VBLLSpec:
    """Architecture + normalisation + provenance of one ported checkpoint (the content of surrogate.json)."""
    d_embed: int = 8
    hidden: int = 64
    n_layers: int = 3
    head_type: str = "het"                 # 'het' (HetRegression) or 'standard' (Regression)
    parameterization: str = "diagonal"
    wishart_scale: float = 1.0
    prior_scale: float = 1.0
    dof: float = 1.0
    noise_prior_scale: float = 0.01
    n_train_per_particle: int = 1         # sets regularization_weight = 1/n; irrelevant at inference
    particles: tuple = PARTICLES
    components: tuple = COMPONENTS
    units: str = "MeV"
    frame: str | None = None              # frame of the training CSV's 4-vectors (None = not yet pinned)
    normaliser: dict = field(default_factory=dict)   # {"input": {particle: {comp: [mean, std]}}, "output": {...}}
    meta: dict = field(default_factory=dict)         # provenance (checkpoint sha256, CSV fingerprint, ...)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["particles"], d["components"] = list(self.particles), list(self.components)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "VBLLSpec":
        d = dict(d)
        d["particles"] = tuple(d.get("particles", PARTICLES)); d["components"] = tuple(d.get("components", COMPONENTS))
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


def build_module(spec: VBLLSpec):
    """The torch module with exactly the trained model's submodule names (so a state_dict loads strictly)."""
    import torch
    import torch.nn as nn
    import vbll

    class ParticleSurrogate(nn.Module):
        def __init__(self):
            super().__init__()
            self.head_type = spec.head_type
            self.embedding = nn.Embedding(2, spec.d_embed)
            d_in = spec.d_embed + 4
            layers = []
            for i in range(spec.n_layers):
                layers += [nn.Linear(d_in if i == 0 else spec.hidden, spec.hidden), nn.ReLU()]
            self.backbone = nn.Sequential(*layers)
            reg_weight = 1.0 / float(spec.n_train_per_particle)
            if spec.head_type == "standard":
                self.reg_heads = nn.ModuleDict({p: vbll.Regression(
                    in_features=spec.hidden, out_features=4, regularization_weight=reg_weight,
                    parameterization=spec.parameterization, prior_scale=spec.prior_scale,
                    wishart_scale=spec.wishart_scale, dof=spec.dof) for p in spec.particles})
            elif spec.head_type == "het":
                # the training code used a subclass (PatchedHetRegression) that only overrides the training
                # loss; parameters and the predictive are those of vbll.HetRegression
                self.reg_heads = nn.ModuleDict({p: vbll.HetRegression(
                    in_features=spec.hidden, out_features=4, regularization_weight=reg_weight,
                    parameterization=spec.parameterization, prior_scale=spec.prior_scale,
                    noise_prior_scale=spec.noise_prior_scale) for p in spec.particles})
            else:
                raise ValueError(f"unknown head_type {spec.head_type!r}")
            self.cls_head = vbll.DiscClassification(spec.hidden, 2, regularization_weight=reg_weight)

        def features(self, type_idx, x_norm):
            return self.backbone(torch.cat([self.embedding(type_idx), x_norm], dim=-1))

        def forward(self, type_idx, truth_4vec):
            """Exactly the trained model's forward (VBLL_SurrogateModel code/model.py): the classification
            head is evaluated first, then one regression head per particle present in the batch, in
            `particles` order — so with the same torch seed the random draws are identical to the
            reference implementation's, which is what the port parity check relies on."""
            h = self.features(type_idx, truth_4vec)
            cls_out = self.cls_head(h)
            reg = {}
            for pidx, pname in enumerate(spec.particles):
                mask = (type_idx == pidx)
                if mask.any():
                    reg[pname] = (mask, self.reg_heads[pname](h[mask]))
            return cls_out, reg

    return ParticleSurrogate()


class VBLLSurrogateModel:
    """A ported checkpoint: spec + torch module. Inputs / outputs in MeV (E, px, py, pz)."""

    def __init__(self, spec: VBLLSpec, module):
        self.spec = spec
        self.module = module.eval()

    # ---- persistence: arrays.npz (state_dict as float32 arrays) + surrogate.json (spec) ------------
    @classmethod
    def from_dir(cls, directory: str | Path) -> "VBLLSurrogateModel":
        import torch
        d = Path(directory)
        js = json.loads((d / "surrogate.json").read_text())
        if js.get("kind") != KIND:
            raise ValueError(f"{d}: kind {js.get('kind')!r} is not {KIND!r}")
        spec = VBLLSpec.from_dict(js["spec"])
        z = np.load(d / "arrays.npz", allow_pickle=False)
        state = {k: torch.from_numpy(np.array(z[k])) for k in z.files}
        module = build_module(spec)
        module.load_state_dict(state, strict=True)
        return cls(spec, module)

    @staticmethod
    def save_dir(directory: str | Path, spec: VBLLSpec, state_dict: dict, extra_json: dict | None = None) -> Path:
        """Write arrays.npz (every tensor of the state_dict as a numpy array) and surrogate.json."""
        d = Path(directory); d.mkdir(parents=True, exist_ok=True)
        arrays = {k: (v.detach().cpu().numpy() if hasattr(v, "detach") else np.asarray(v)) for k, v in state_dict.items()}
        np.savez_compressed(d / "arrays.npz", **arrays)
        js = {"kind": KIND, "binning": None, "spec": spec.to_dict(),
              "arrays": {k: [list(v.shape), str(v.dtype)] for k, v in arrays.items()}}
        js.update(extra_json or {})
        (d / "surrogate.json").write_text(json.dumps(js, indent=2, default=str))
        return d

    # ---- normalisation ------------------------------------------------------------------------------
    def _stats(self, kind: str, particle: str) -> tuple[np.ndarray, np.ndarray]:
        tab = self.spec.normaliser[kind][particle]
        mu = np.array([tab[c][0] for c in self.spec.components], dtype=np.float32)
        sd = np.array([tab[c][1] for c in self.spec.components], dtype=np.float32)
        return mu, sd

    def normalise_input(self, particle: str, truth_mev: np.ndarray) -> np.ndarray:
        mu, sd = self._stats("input", particle)
        return (np.asarray(truth_mev, np.float32) - mu) / (sd + np.float32(1e-8))

    def denormalise_output(self, particle: str, v_norm: np.ndarray) -> np.ndarray:
        mu, sd = self._stats("output", particle)
        return np.asarray(v_norm, np.float32) * (sd + np.float32(1e-8)) + mu

    def denormalise_input(self, particle: str, v_norm: np.ndarray) -> np.ndarray:
        mu, sd = self._stats("input", particle)
        return np.asarray(v_norm, np.float32) * (sd + np.float32(1e-8)) + mu

    def output_scale(self, particle: str) -> np.ndarray:
        return self._stats("output", particle)[1] + np.float32(1e-8)

    # ---- inference ----------------------------------------------------------------------------------
    def predictive_normalised(self, particle: str, x_norm, generator=None):
        """One draw of the head's predictive Normal (mean, sigma) in normalised units, exactly as the
        reference `evaluate.collect_predictions` obtains it (`reg_out.predictive`)."""
        import torch
        x = torch.as_tensor(np.asarray(x_norm, np.float32))
        t = torch.full((x.shape[0],), PARTICLE_INDEX[particle], dtype=torch.long)
        with torch.no_grad():
            h = self.module.features(t, x)
            head = self.module.reg_heads[particle]
            if self.spec.head_type == "het":
                if generator is not None:
                    state = torch.random.get_rng_state(); torch.random.set_rng_state(generator.get_state())
                dist = head.predictive_sample(h, consistent_variance=False)
                if generator is not None:
                    generator.set_state(torch.random.get_rng_state()); torch.random.set_rng_state(state)
            else:
                dist = head.predictive(h)
        return dist.mean.numpy(), dist.stddev.numpy()

    def predict(self, particle: str, truth_mev: np.ndarray, generator=None) -> tuple[np.ndarray, np.ndarray]:
        """(mean, sigma) of the reconstructed 4-vector in MeV for truth 4-vectors in MeV, shape (n, 4)."""
        m, s = self.predictive_normalised(particle, self.normalise_input(particle, truth_mev), generator)
        return self.denormalise_output(particle, m), s * self.output_scale(particle)

    def smear(self, particle: str, truth_mev: np.ndarray, n_samples: int = 1, seed: int = 0,
              batch_size: int = 65536) -> np.ndarray:
        """Reconstructed 4-vectors in MeV, shape (n_samples, n, 4): per sample a fresh sigma draw
        (het head) and a Normal(mean, sigma) draw — the model's own predictive marginal."""
        import torch
        g = torch.Generator().manual_seed(int(seed))
        truth_mev = np.asarray(truth_mev, np.float32)
        out = np.empty((n_samples, len(truth_mev), 4), dtype=np.float32)
        for k in range(n_samples):
            for lo in range(0, len(truth_mev), batch_size):
                m, s = self.predict(particle, truth_mev[lo:lo + batch_size], generator=g)
                eps = torch.randn(m.shape, generator=g).numpy()
                out[k, lo:lo + batch_size] = m + s * eps
        return out


def load_vbll_model(directory: str | Path) -> VBLLSurrogateModel:
    return VBLLSurrogateModel.from_dir(directory)
