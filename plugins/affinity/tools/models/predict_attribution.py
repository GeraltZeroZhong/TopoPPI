#!/usr/bin/env python3
"""Apply one saved attribution model to a label-free named feature matrix."""
import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from topoppi_affinity.benchmark_models import LoadedFeatures
from topoppi_affinity.thermodynamics import pkd_to_dg, pkd_to_kd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--features', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    bundle = joblib.load(args.model)
    features = LoadedFeatures.read(args.features)
    lookup = {name: i for i, name in enumerate(features.feature_names)}
    x = features.X[:, [lookup[name] for name in bundle['feature_names']]]
    usable = np.isfinite(x).all(axis=1)
    values = np.full(len(x), np.nan)
    if usable.any():
        values[usable] = bundle['model'].predict(x[usable])
    usable &= np.isfinite(values)
    result = pd.DataFrame(dict(sample_id=features.sample_ids, predicted_pkd=values,
                               status=np.where(usable, 'ok', 'error')))
    result['predicted_kd_M'], result['predicted_dg_298_15_kcal_mol'] = np.nan, np.nan
    result.loc[usable, 'predicted_kd_M'] = pkd_to_kd(values[usable])
    result.loc[usable, 'predicted_dg_298_15_kcal_mol'] = pkd_to_dg(values[usable], 298.15)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f'{int(usable.sum())}/{len(result)} predictions; no labels read.')
