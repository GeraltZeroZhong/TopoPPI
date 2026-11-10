"""Load portable linear and nonlinear models for frozen inference."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .nonlinear_ranking import TorchPairwiseMLP, XGBPairwiseRanker
from .representations import center_within_queries


def predict_frozen(model_dir, frame, embeddings=None):
    """No labels, fitting or model selection occur in this function."""
    model_dir = Path(model_dir)
    metadata = json.loads((model_dir / 'model.json').read_text())
    X = frame[metadata['columns']].to_numpy(float)
    if (model_dir / 'embedding_pca.npz').exists():
        if embeddings is None:
            raise ValueError('This frozen model requires independent WT ESMC embeddings')
        centered = center_within_queries(embeddings, frame.query_id.to_numpy())
        with np.load(model_dir / 'embedding_pca.npz', allow_pickle=False) as pca:
            projected = (centered - pca['mean']) @ pca['components'].T
        X = np.column_stack([X, projected])
    if metadata['kind'] == 'linear':
        with np.load(model_dir / 'linear_model.npz', allow_pickle=False) as model:
            if "median" in model:
                missing = np.isnan(X)
                X = np.where(missing, model["median"], X)
                indicators = model["missing_indicator_indices"].astype(int)
                if len(indicators):
                    X = np.column_stack([X, missing[:, indicators]])
            scores = ((X - model['mean']) / model['scale']) @ model['coefficient']
        return scores
    cls = {'xgb': XGBPairwiseRanker, 'mlp': TorchPairwiseMLP}[metadata['kind']]
    return np.mean([cls.load(model_dir / f'seed_{seed}').decision_function(X) for seed in metadata['seeds']], axis=0)
