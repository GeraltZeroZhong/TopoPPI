"""Optional mass-weighted finite graph readouts with validation stopping."""

import copy
from dataclasses import asdict
from time import perf_counter

import numpy as np
import torch
from torch import nn

from .graphs import FEATURE_NAMES, TARGET_NAMES, ReaderConfig


def _batch(samples, view, device="cpu"):
    features, mass, pool_index, edges = [], [], [], []
    offset = 0
    for batch_id, sample in enumerate(samples):
        if view not in sample.edges:
            raise ValueError(f"Missing {view} graph for {sample.sample_id}; do not silently omit it")
        features.append(sample.features)
        mass.append(sample.mass)
        pool_index.append(2 * batch_id + sample.side)
        edges.append(sample.edges[view] + offset)
        offset += len(sample.features)
    pairs = np.concatenate(edges)
    directed = np.concatenate([pairs, pairs[:, ::-1]])
    result = {
        "x": torch.as_tensor(np.concatenate(features), dtype=torch.float32),
        "mass": torch.as_tensor(np.concatenate(mass), dtype=torch.float32),
        "pool_index": torch.as_tensor(np.concatenate(pool_index), dtype=torch.long),
        "edges": torch.as_tensor(directed, dtype=torch.long),
        "size": len(samples),
        "target": torch.as_tensor(np.array([s.target for s in samples]), dtype=torch.float32),
        "mask": torch.as_tensor(np.array([s.target_mask for s in samples]), dtype=torch.bool),
    }
    return {k: v.to(device) if isinstance(v, torch.Tensor) else v for k, v in result.items()}


INPUT_CHANNELS = len(FEATURE_NAMES)
DEFAULT_READER_CONFIG = ReaderConfig()


class FiniteGraphReader(nn.Module):
    """The same mass-weighted message operator and two-side pooling in all arms."""

    def __init__(self, width=32, depth=1, input_channels=INPUT_CHANNELS):
        super().__init__()
        self.input_layer = nn.Linear(input_channels, width)
        self.self_layers = nn.ModuleList([nn.Linear(width, width) for _ in range(depth)])
        self.neighbor_layers = nn.ModuleList([nn.Linear(width, width, bias=False) for _ in range(depth)])
        self.o_head = nn.Sequential(nn.Linear(2 * width, width), nn.ReLU(), nn.Linear(width, 1))
        self.x_head = nn.Sequential(nn.Linear(2 * width, width), nn.ReLU(), nn.Linear(width, 11))

    def forward(self, batch):
        h = torch.relu(self.input_layer(batch["x"]))
        source, destination = batch["edges"].T
        for own, neighbors in zip(self.self_layers, self.neighbor_layers, strict=False):
            summed = torch.zeros_like(h)
            norm = torch.zeros((len(h), 1), dtype=h.dtype, device=h.device)
            if len(source):
                weights = batch["mass"][source, None]
                summed.index_add_(0, destination, h[source] * weights)
                norm.index_add_(0, destination, weights)
            h = torch.relu(own(h) + neighbors(summed / norm.clamp_min(1e-12)))
        pooled = torch.zeros((2 * batch["size"], h.shape[1]), dtype=h.dtype, device=h.device)
        norm = torch.zeros((2 * batch["size"], 1), dtype=h.dtype, device=h.device)
        pooled.index_add_(0, batch["pool_index"], h * batch["mass"][:, None])
        norm.index_add_(0, batch["pool_index"], batch["mass"][:, None])
        pooled = (pooled / norm.clamp_min(1e-12)).reshape(batch["size"], 2, -1)
        focal_pair = torch.cat([pooled, pooled.flip(dims=(1,))], dim=-1)
        o = torch.sigmoid(self.o_head(focal_pair).squeeze(-1))
        pair = torch.cat([pooled.sum(dim=1), torch.abs(pooled[:, 0] - pooled[:, 1])], dim=-1)
        x = torch.softmax(self.x_head(pair), dim=-1)
        return torch.cat([o, x], dim=-1)


def reader_loss(prediction, target, mask):
    """Equal O/X-family loss: O MSE and half squared-L2 distribution error.

    O sides are averaged within each sample; X is one complete distribution.
    Both targets are bounded, so no train/test-fitted target scaling is needed.
    """
    target = torch.where(mask, target, torch.zeros_like(target))
    o_mask = mask[:, :2]
    o_count = o_mask.sum(dim=1)
    o_error = (((prediction[:, :2] - target[:, :2]) ** 2) * o_mask).sum(dim=1) / o_count.clamp_min(1)
    x_valid = mask[:, 2:].all(dim=1)
    x_error = ((prediction[:, 2:] - target[:, 2:]) ** 2).sum(dim=1) / 2
    count = (o_count > 0).float() + x_valid.float()
    error = (o_error * (o_count > 0) + x_error * x_valid) / count.clamp_min(1)
    valid = count > 0
    return error[valid].mean() if torch.any(valid) else prediction.sum() * 0


def predict_reader(model, samples, view, batch_size=8):
    model.eval()
    outputs = []
    with torch.no_grad():
        for start in range(0, len(samples), batch_size):
            outputs.append(
                model(_batch(samples[start : start + batch_size], view, next(model.parameters()).device)).cpu().numpy()
            )
    return np.concatenate(outputs) if outputs else np.empty((0, 13), dtype=np.float32)


def load_reader_checkpoint(path):
    payload = torch.load(path, map_location="cpu", weights_only=True)
    if tuple(payload["feature_names"]) != FEATURE_NAMES or tuple(payload["target_names"]) != TARGET_NAMES:
        raise ValueError("Checkpoint feature/target contract differs from current reader")
    config = payload["config"]
    model = FiniteGraphReader(width=config["width"], depth=config["depth"])
    model.load_state_dict(payload.pop("state_dict"))
    model.eval()
    return model, payload


def fit_reader(training_samples, validation_samples, view, config=DEFAULT_READER_CONFIG):
    """Fit with validation-only stopping; outer test data is not accepted."""
    if not training_samples or not validation_samples:
        raise ValueError("Separate nonempty training and inner-validation groups are required")
    train_groups = {s.family_id for s in training_samples}
    if train_groups & {s.family_id for s in validation_samples}:
        raise ValueError("Training and validation relationship groups overlap")
    started = perf_counter()
    torch.set_num_threads(config.cpu_threads)
    torch.manual_seed(config.seed)
    if config.device.startswith("cuda"):
        torch.cuda.manual_seed_all(config.seed)
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.cuda.reset_peak_memory_stats()
    rng = np.random.default_rng(config.seed)
    model = FiniteGraphReader(config.width, config.depth).to(config.device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay)
    best, best_state, best_epoch, stale, history = float("inf"), None, 0, 0, []
    validation = _batch(validation_samples, view, config.device)
    for epoch in range(config.max_epochs):
        model.train()
        order = rng.permutation(len(training_samples))
        losses = []
        for start in range(0, len(order), config.batch_size):
            batch = _batch([training_samples[i] for i in order[start : start + config.batch_size]], view, config.device)
            optimizer.zero_grad()
            loss = reader_loss(model(batch), batch["target"], batch["mask"])
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach()))
        model.eval()
        with torch.no_grad():
            validation_loss = float(reader_loss(model(validation), validation["target"], validation["mask"]))
        history.append({"epoch": epoch + 1, "train_loss": float(np.mean(losses)), "validation_loss": validation_loss})
        if validation_loss < best:
            best, best_epoch, stale = validation_loss, epoch + 1, 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1
        if stale >= config.patience:
            break
    if best_state is None:
        raise RuntimeError("No finite validation checkpoint was produced")
    model.load_state_dict(best_state)
    if config.device.startswith("cuda"):
        torch.cuda.synchronize()
    return model, {
        "config": asdict(config),
        "view": view,
        "best_epoch": best_epoch,
        "best_validation_loss": best,
        "history": history,
        "fit_seconds": perf_counter() - started,
        "nparams": sum(p.numel() for p in model.parameters()),
        "parameter_bytes": sum(p.numel() * p.element_size() for p in model.parameters()),
        "train_groups": sorted(train_groups),
        "validation_groups": sorted({s.family_id for s in validation_samples}),
        "target_scaling": "none: O and X are bounded [0,1]",
        "device": config.device,
        "peak_cuda_bytes": torch.cuda.max_memory_allocated() if config.device.startswith("cuda") else None,
        "cuda_scatter_note": "Float32 atomic reductions may vary at rounding precision across executions"
        if config.device.startswith("cuda")
        else None,
        "hyperparameter_candidates": 1,
    }
