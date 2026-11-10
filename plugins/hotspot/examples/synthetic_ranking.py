"""Fit and score two invented residue queries using an analytical toy feature.

The eight feature rows and labels are constructed here without structural or
experimental records. Each query has two positives and two negatives.
"""
import numpy as np
from topoppi_hotspot.ranking_metrics import ranking_metrics
from topoppi_hotspot.ranking_models import PairwiseHotspotRanker

features = np.array([[-2.], [-1.], [1.], [2.], [-3.], [-1.], [1.], [3.]])
labels = np.array([0, 0, 1, 1, 0, 0, 1, 1])
queries = np.array(["toy_a"] * 4 + ["toy_b"] * 4)
model = PairwiseHotspotRanker(alpha=0.1).fit(features, labels, queries)
scores = model.decision_function(features)
for query in np.unique(queries):
    selected = queries == query
    values = ranking_metrics(labels[selected], scores[selected])
    print(f"{query}: scores={scores[selected].round(3).tolist()}, AP={values['AP']:.1f}")
