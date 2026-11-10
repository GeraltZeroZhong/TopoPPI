"""Weighted kernel prediction for bilinear features."""
from dataclasses import dataclass


@dataclass
class WeightedKernelRegressor:
    estimator: object
    blocks: tuple
    weights: tuple
    center: float

    def predict(self, matrix):
        gram = sum(weight * block.transform(matrix) for block, weight in
                   zip(self.blocks, self.weights, strict=True) if weight > 0)
        return self.estimator.predict(gram) + self.center
