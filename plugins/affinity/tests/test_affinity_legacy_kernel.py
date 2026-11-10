import numpy as np

from topoppi_affinity.legacy_kernel import concatenate_blocks


def test_contact_ablation_removes_surface_from_both_blocks():
    plm = np.arange(6).reshape(2, 3) + 1.
    legacy = np.arange(8).reshape(2, 4) + 20.
    pnames = ['plm.global', 'surface.weighted', 'plm.contact']
    lnames = ['x', 'hidden_surface_name', 'y', 'z']
    groups = ['contact', 'surface', 'sequence', 'contact']
    full, names, second = concatenate_blocks(plm, pnames, legacy, lnames, groups, use_surface=True)
    control, _, control_second = concatenate_blocks(plm, pnames, legacy, lnames, groups, use_surface=False)
    np.testing.assert_array_equal(control[:, [1, 4]], 0.)
    np.testing.assert_array_equal(control[:, [0, 2, 3, 5, 6]], full[:, [0, 2, 3, 5, 6]])
    np.testing.assert_array_equal(second, control_second)
    np.testing.assert_array_equal(full[:, ~second], plm)
    np.testing.assert_array_equal(full[:, second], legacy)
    assert len(set(names)) == 7
    # Input caches remain unchanged for the subsequent surface run.
    assert plm[0, 1] == 2. and legacy[0, 1] == 21.


def test_bundle_roundtrip_unlabeled_reordered_ids_and_units(tmp_path):
    import joblib
    import pandas as pd
    import pytest

    from topoppi_affinity.legacy_kernel import fit_legacy_bundle, predict_legacy_kernel
    from topoppi_affinity.thermodynamics import R_KCAL

    generator = np.random.default_rng(42)
    plm, legacy = generator.normal(size=(8, 3)), generator.normal(size=(8, 4))
    pnames = np.array(['plm.a', 'surface.a', 'plm.b'])
    lnames = np.array(['chem.a', 'area', 'chem.b', 'count'])
    groups = np.array(['contact', 'surface', 'sequence', 'contact'])
    frame = pd.DataFrame({'sample_id': [f's{i}' for i in range(8)], 'complex_id': [f'c{i}' for i in range(8)],
                          'family_group': [f'f{i//2}' for i in range(8)], 'pkd': generator.normal(7., 1., 8)})
    bundle = fit_legacy_bundle(plm, pnames, legacy, lnames, groups, frame)
    path = tmp_path / 'model.joblib'
    joblib.dump(bundle, path)
    first, second = tmp_path / 'plm.npz', tmp_path / 'legacy.npz'
    np.savez(first, sample_ids=frame.sample_id.to_numpy(str), X=plm, feature_names=pnames)
    np.savez(second, sample_ids=frame.sample_id.to_numpy(str)[::-1], X=legacy[::-1], feature_names=lnames)
    result = predict_legacy_kernel(path, first, second, temperature_K=300.)
    combined, _, _ = concatenate_blocks(plm, pnames, legacy, lnames, groups, use_surface=True)
    np.testing.assert_allclose(result.predicted_pkd, bundle['model'].predict(combined), rtol=0, atol=1e-12)
    assert result.sample_id.tolist() == frame.sample_id.tolist()
    np.testing.assert_allclose(result.predicted_kd_M, 10.**(-result.predicted_pkd))
    np.testing.assert_allclose(result.predicted_dg_at_reference_kcal_mol, -R_KCAL*300.*np.log(10)*result.predicted_pkd)
    np.savez(second, sample_ids=frame.sample_id.to_numpy(str), X=legacy, feature_names=lnames[::-1])
    with pytest.raises(ValueError, match='columns differ'):
        predict_legacy_kernel(path, first, second)
