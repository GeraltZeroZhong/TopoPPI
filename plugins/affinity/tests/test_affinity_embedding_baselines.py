"""Scientific invariants and complete surface ablation in frozen representations."""
import joblib
import numpy as np
import pytest
import torch
from torch_geometric.data import Data

from topoppi_affinity.embedding_baselines import fit_frozen_regressor, pooled_features, zero_surface_slots


@pytest.fixture()
def paired():
    generator = torch.Generator().manual_seed(13)
    x = torch.rand((6, 45), generator=generator)
    x_mut = x.clone()
    x_mut[[1, 4], :26] += .5
    esm = torch.rand((6, 4), generator=generator)
    esm_mut = esm.clone()
    esm_mut[[1, 4]] += .2
    return Data(x=x, x_mut=x_mut, esm=esm, esm_mut=esm_mut,
                side=torch.tensor([0, 0, 0, 1, 1, 1]), contact_mask=torch.tensor([True, True, False, True, True, False]),
                mutation_mask=torch.tensor([False, True, False, False, True, False]))


@pytest.mark.parametrize("surface", [False, True])
def test_partner_symmetry_and_joint_mutation_oddness(paired, surface):
    absolute, _ = pooled_features(paired, target="pkd", use_surface=surface)
    exchanged = paired.clone()
    exchanged.side = 1 - paired.side
    np.testing.assert_allclose(absolute, pooled_features(exchanged, target="pkd", use_surface=surface)[0])
    mutation, _ = pooled_features(paired, target="ddg", use_surface=surface)
    reverse = paired.clone()
    reverse.x, reverse.x_mut = paired.x_mut.clone(), paired.x.clone()
    reverse.esm, reverse.esm_mut = paired.esm_mut.clone(), paired.esm.clone()
    np.testing.assert_allclose(mutation, -pooled_features(reverse, target="ddg", use_surface=surface)[0])
    noop = paired.clone()
    noop.x_mut, noop.esm_mut = noop.x.clone(), noop.esm.clone()
    assert not pooled_features(noop, target="ddg", use_surface=surface)[0].any()


@pytest.mark.parametrize("target", ["pkd", "ddg"])
def test_no_surface_pooling_ignores_all_surface_columns(paired, target):
    original, names = pooled_features(paired, target=target, use_surface=False)
    paired.x[:, 30:] += 1000.
    paired.x_mut[:, 30:] += 1000.
    np.testing.assert_array_equal(original, pooled_features(paired, target=target, use_surface=False)[0])
    assert not any(name.startswith("surface.") for name in names)


@pytest.mark.parametrize("family", ["ridge", "rbf", "catboost"])
def test_fitted_odd_regressor_and_saved_model(family, tmp_path):
    x = np.array([[1., 2.], [2., -1.], [3., 1.], [-2., -1.]])
    y = np.array([1., 2., 4., -1.])
    params = {"alpha": .1, "gamma": .5} if family != "catboost" else {"depth": 2, "iterations": 20, "thread_count": 8}
    model = fit_frozen_regressor(x, y, ["a", "a", "b", "c"], family=family, params=params, antisymmetric=True)
    if family == "catboost":
        assert model.estimator.get_params()['thread_count'] == 8
    np.testing.assert_allclose(model.predict(x), -model.predict(-x), atol=1e-12)
    np.testing.assert_array_equal(model.predict(np.zeros((1, 2))), [0.])
    if model.scaler is not None:
        np.testing.assert_allclose(model.scaler.mean_, 0., atol=1e-12)
    path = tmp_path / "model.joblib"
    joblib.dump(model, path)
    np.testing.assert_array_equal(joblib.load(path).predict(x), model.predict(x))


@pytest.mark.parametrize("target", ["pkd", "ddg"])
def test_contact_control_preserves_all_non_surface_entries_and_dimension(paired, target):
    full, names = pooled_features(paired, target=target, use_surface=True)
    control = zero_surface_slots(full, names)
    original, _ = pooled_features(paired, target=target, use_surface=False)
    surface = np.char.startswith(names, "surface.")
    assert full.shape == control.shape
    np.testing.assert_array_equal(control[~surface], full[~surface])
    np.testing.assert_array_equal(control[~surface], original)
    assert not control[surface].any()


def test_normalized_kernel_matches_training_transform_and_is_psd():
    from topoppi_affinity.embedding_baselines import fit_kernel_block

    values = np.random.default_rng(7).normal(size=(12, 6))
    columns = np.array([True, True, True, False, False, False])
    block, gram = fit_kernel_block(values, columns, np.ones(12), .1)
    np.testing.assert_allclose(gram, block.transform(values), atol=1e-12)
    np.testing.assert_allclose(gram.mean(axis=0), 0., atol=1e-12)
    assert np.isclose(np.trace(gram) / len(gram), 1.)
    assert np.linalg.eigvalsh(gram).min() > -1e-10


def test_multikernel_oddness_and_zero_surface_weight_control(tmp_path):
    from sklearn.kernel_ridge import KernelRidge

    from topoppi_affinity.embedding_baselines import MultiKernelRegressor, fit_kernel_block

    original = np.random.default_rng(31).normal(size=(8, 6))
    train = np.concatenate([original, -original])
    labels = np.concatenate([np.arange(8), -np.arange(8)])
    surface = np.array([False, False, False, True, True, True])
    plm_block, gram = fit_kernel_block(train, ~surface, np.ones(16), .1)
    surface_block, _ = fit_kernel_block(train, surface, np.ones(16), .1)
    estimator = KernelRidge(alpha=.1, kernel="precomputed").fit(gram, labels)
    model = MultiKernelRegressor(estimator, (plm_block, surface_block), 0., True, 0.)
    prediction = model.predict(original)
    np.testing.assert_allclose(prediction, -model.predict(-original), atol=1e-12)
    np.testing.assert_array_equal(model.predict(np.zeros((1, 6))), [0.])
    changed = original.copy()
    changed[:, surface] += 1000.
    np.testing.assert_array_equal(prediction, model.predict(changed))
    path = tmp_path / "multikernel.joblib"
    joblib.dump(model, path)
    np.testing.assert_array_equal(prediction, joblib.load(path).predict(original))
