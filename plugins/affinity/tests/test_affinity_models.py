import numpy as np
import pytest

from topoppi_affinity.evaluation import regression_metrics
from topoppi_affinity.models import complex_weights, fit_regressor
from topoppi_affinity.thermodynamics import ddg_from_kd, kd_to_pkd, pkd_to_dg, pkd_to_kd


def test_thermodynamic_units_sign_and_roundtrip():
    assert kd_to_pkd(1e-9) == pytest.approx(9)
    assert pkd_to_kd(9) == pytest.approx(1e-9)
    assert pkd_to_dg(9) == pytest.approx(-12.278223, abs=1e-5)
    assert ddg_from_kd(1e-9, 1e-8, 298.15) == pytest.approx(1.364247, abs=1e-5)
    assert ddg_from_kd(1e-8, 1e-9, 298.15) < 0
    with pytest.raises(ValueError):
        kd_to_pkd(0)
    with pytest.raises(ValueError):
        pkd_to_dg(8, -1)


def test_equal_complex_weight_and_macro_metrics():
    groups = np.array(["a", "a", "a", "b"])
    weights = complex_weights(groups)
    assert weights[:3].sum() == pytest.approx(weights[3])
    report = regression_metrics([0, 0, 0, 0], [1, 1, 1, 3], groups, n_bootstrap=30)
    assert report["mae"] == pytest.approx(1.5)
    assert report["macro_mae"] == pytest.approx(2)
    assert report["pearson"] is None
    assert report == regression_metrics([0, 0, 0, 0], [1, 1, 1, 3], groups, n_bootstrap=30)


def test_training_only_scaler_schema_and_saved_predictions(tmp_path):
    import joblib

    matrix = np.array([[0., 1.], [1., 2.], [2., 3.], [3., 4.]])
    model = fit_regressor(matrix, [5., 6., 7., 8.], ["a", "b"], target="pkd",
                          model_name="ridge_1", groups=["1", "2", "3", "4"])
    assert model.estimator.named_steps["scale"].mean_ == pytest.approx([1.5, 2.5])
    predicted = model.predict([[100., 101.]], ["a", "b"])
    assert model.estimator.named_steps["scale"].mean_ == pytest.approx([1.5, 2.5])
    assert model.predict([[101., 100.]], ["b", "a"]) == pytest.approx(predicted)
    with pytest.raises(ValueError, match="Missing"):
        model.predict([[100.]], ["a"])
    model.save(tmp_path / "model.joblib")
    restored = joblib.load(tmp_path / "model.joblib")
    assert restored.predict([[100., 101.]], ["a", "b"]) == pytest.approx(predicted)


def test_noop_mutation_has_exact_zero_prediction():
    model = fit_regressor([[1.], [2.], [3.]], [2., 3., 4.], ["mutation"], target="ddg",
                          model_name="ridge_1", groups=["a", "a", "b"])
    assert model.predict([[0.]], ["mutation"])[0] == 0.


def test_noop_zero_with_nonzero_wild_type_context():
    model = fit_regressor([[8., 1.], [9., 1.], [10., 2.]], [2., 3., 4.],
                          ["contact.area", "mutation.n_mutations"], target="ddg",
                          model_name="ridge_1", groups=["a", "b", "c"])
    assert model.predict([[10., 0.]], ["contact.area", "mutation.n_mutations"])[0] == 0.
