"""Existing GINE/GATv2 and APPNP operators adapted to TopoPPI residue graphs.

GINE: Hu et al., ICLR 2020; GATv2: Brody et al., ICLR 2022;
APPNP: Klicpera et al., ICLR 2019. Implementations are PyTorch Geometric's
original operators. This is an input/readout adaptation, not a novel algorithm.
WT/MT share the encoder. An antisymmetric pair readout models local changes;
absolute affinity uses symmetric two-partner pooling.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import torch
from torch import nn
from torch_geometric.nn import APPNP, GATv2Conv, GINEConv, global_add_pool


@dataclass
class GraphModelConfig:
    architecture: str = 'gine'
    hidden: int = 96
    layers: int = 3
    dropout: float = .15
    use_plm: bool = True
    use_surface: bool = True
    target: str = 'pkd'
    plm_dimension: int = 1152
    node_dimension: int = 45
    max_nodes: int | None = 512
    diffusion_root: str = 'diffusion'
    surface_width: int = 32
    surface_blocks: int = 2

    def to_dict(self):
        return asdict(self)


def physical_transform(x):
    return torch.cat([x[:, :21], torch.sign(x[:, 21:]) * torch.log1p(x[:, 21:].abs())], dim=-1)


def edge_transform(edge_attr):
    return torch.cat([edge_attr[:, :2] / 10., edge_attr[:, 2:5],
                      torch.log1p(edge_attr[:, 5:].clamp_min(0.)) / 3.], dim=-1)


@lru_cache(maxsize=4)
def _surface_operators(path):
    from topoppi_affinity.surface_diffusion import load_surface_operators
    return load_surface_operators(path)


class AffinityGraphNetwork(nn.Module):
    def __init__(self, config: GraphModelConfig):
        super().__init__()
        self.config = config
        d = config.hidden
        self.register_buffer('feature_center', torch.zeros(config.node_dimension))
        self.register_buffer('feature_scale', torch.ones(config.node_dimension))
        self.register_buffer('target_center', torch.tensor(0.))
        self.register_buffer('target_scale', torch.tensor(1.))
        self.chemistry = nn.Linear(config.node_dimension, d)
        self.language = nn.Sequential(nn.LayerNorm(config.plm_dimension), nn.Linear(config.plm_dimension, d))
        self.input_norm = nn.LayerNorm(d)
        self.layers = nn.ModuleList()
        self.norms = nn.ModuleList([nn.LayerNorm(d) for _ in range(config.layers)])
        for _ in range(config.layers):
            if config.architecture in {'gine', 'diffusion_gine'}:
                self.layers.append(GINEConv(nn.Sequential(nn.Linear(d, d), nn.SiLU(), nn.Linear(d, d)),
                                           edge_dim=8, train_eps=True))
            elif config.architecture == 'gatv2':
                self.layers.append(GATv2Conv(d, d // 4, heads=4, edge_dim=8, dropout=config.dropout))
            elif config.architecture in {'appnp', 'mlp'}:
                self.layers.append(nn.Sequential(nn.Linear(d, d), nn.SiLU(), nn.Linear(d, d)))
            else:
                raise ValueError(f'Unknown existing operator: {config.architecture}')
        self.diffusion = APPNP(K=10, alpha=.1, dropout=0.)
        if config.architecture == 'diffusion_gine':
            from topoppi_affinity.surface_diffusion import SurfaceEncoder
            self.surface_encoder = SurfaceEncoder(d, d, width=config.surface_width, n_blocks=config.surface_blocks)
            self.surface_fusion = nn.Linear(2 * d, d)

        self.dropout = nn.Dropout(config.dropout)
        self.absolute_head = nn.Sequential(nn.Linear(2 * d + 2, d), nn.SiLU(), nn.Dropout(config.dropout), nn.Linear(d, 1))
        self.mutation_head = nn.Sequential(nn.Linear(3 * d + 1, d), nn.SiLU(), nn.Linear(d, d // 2),
                                           nn.SiLU(), nn.Linear(d // 2, 1))

    def encode(self, data, mutant=False):
        x = data.x_mut if mutant else data.x
        x = (physical_transform(x) - self.feature_center) / self.feature_scale
        h = self.chemistry(x)
        if self.config.use_plm:
            embedding = data.esm_mut if mutant else data.esm
            h = h + self.language(embedding.float())
        h = self.input_norm(h)
        if self.config.architecture == 'diffusion_gine':
            from topoppi_affinity.surface_diffusion import operators_to
            if data.num_graphs != 1:
                raise ValueError('Surface diffusion uses full graphs, microbatch size one and gradient accumulation')
            path = Path(self.config.diffusion_root) / data.source[0] / (data.complex_id[0].replace(':', '_').replace(',', '+') + '.pt')
            operators = operators_to(_surface_operators(str(path)), h.device)
            surface = self.surface_encoder(h, operators)
            h = h + self.surface_fusion(torch.cat([surface['surface'], surface['interface']], dim=-1))
        edges = edge_transform(data.edge_attr)
        for layer, norm in zip(self.layers, self.norms, strict=True):
            update = layer(h) if self.config.architecture in {'appnp', 'mlp'} else layer(h, data.edge_index, edges)
            h = norm(h + self.dropout(torch.nn.functional.silu(update)))
        if self.config.architecture == 'appnp':
            # Standard personalized PageRank diffusion. Surface boundary weights
            # distinguish actual footprint adjacency from the contact graph.
            weight = torch.exp(-data.edge_attr[:, 0] / 10.)
            if self.config.use_surface:
                weight = weight + torch.log1p(data.edge_attr[:, 6])
            h = self.diffusion(h, data.edge_index, edge_weight=weight)
        return h

    def absolute_readout(self, h, data):
        groups = data.batch * 2 + data.side
        weights = data.contact_mask.float().reshape(-1, 1)
        n_graphs = data.num_graphs
        mass = global_add_pool(weights, groups, size=2 * n_graphs)
        pooled = global_add_pool(h * weights, groups, size=2 * n_graphs) / mass.clamp_min(1.)
        paired = pooled.reshape(n_graphs, 2, -1)
        count = torch.log1p(mass.reshape(n_graphs, 2))
        features = torch.cat([paired.sum(1), (paired[:, 0] - paired[:, 1]).abs(),
                              count.sum(1, keepdim=True), (count[:, :1] - count[:, 1:]).abs()], dim=-1)
        return self.absolute_head(features).reshape(-1)

    def forward(self, data):
        # The same stochastic encoder represents both molecular states. Sharing
        # dropout masks prevents a nonzero training signal for an identity pair.
        paired_dropout = self.training and self.config.target == 'ddg'
        if paired_dropout:
            cpu_rng = torch.get_rng_state()
            devices = [data.x.device.index] if data.x.is_cuda else []
            cuda_rng = torch.cuda.get_rng_state(data.x.device) if devices else None
        wt = self.encode(data)
        if self.config.target == 'pkd':
            return self.absolute_readout(wt, data) * self.target_scale + self.target_center
        if paired_dropout:
            with torch.random.fork_rng(devices=devices):
                torch.set_rng_state(cpu_rng)
                if devices:
                    torch.cuda.set_rng_state(cuda_rng, data.x.device)
                mt = self.encode(data, mutant=True)
        else:
            mt = self.encode(data, mutant=True)
        n_graphs = data.num_graphs
        mask = data.mutation_mask.float().reshape(-1, 1)
        delta = global_add_pool((mt - wt) * mask, data.batch, size=n_graphs)
        count = global_add_pool(mask, data.batch, size=n_graphs)
        context = global_add_pool((mt + wt) * .5 * mask, data.batch, size=n_graphs) / count.clamp_min(1.)
        global_delta = global_add_pool((mt - wt) * data.contact_mask[:, None], data.batch, size=n_graphs)
        interface_count = global_add_pool(data.contact_mask[:, None].float(), data.batch, size=n_graphs)
        global_delta = global_delta / interface_count.clamp_min(1.)
        forward = torch.cat([delta, context, global_delta, torch.log1p(count)], dim=-1)
        reverse = torch.cat([-delta, context, -global_delta, torch.log1p(count)], dim=-1)
        prediction = .5 * (self.mutation_head(forward) - self.mutation_head(reverse)).reshape(-1)
        return prediction * self.target_scale
