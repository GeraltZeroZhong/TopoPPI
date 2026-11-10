"""ID alignment and explicit prior-failure behavior of fixed model fusion."""
import joblib
import numpy as np
import pandas as pd
import pytest

from topoppi_affinity.embedding_baselines import fit_frozen_regressor
from topoppi_affinity.mutation_ensemble import bundle_mutation_ensemble, predict_mutation_ensemble
from topoppi_affinity.residual_baselines import ResidualRegressor


def test_serialized_equal_fusion_and_base_only_fallback(tmp_path):
    x = np.arange(24, dtype=float).reshape(8, 3) / 10.
    ids, names = np.asarray([f's{i}' for i in range(8)]), np.asarray(['plm.a','contact.b','surface.c'])
    base = fit_frozen_regressor(x, x[:,0] * .8, np.repeat(['a','b'],4),
                               family='ridge',params={'alpha':1.})
    residual = ResidualRegressor(1.2, fit_frozen_regressor(x, x[:,1] * .2,np.repeat(['a','b'],4),
                                family='ridge',params={'alpha':1.},antisymmetric=True))
    joblib.dump({'model':base,'feature_names':names,'report':{'target':'ddg'}},tmp_path/'base.joblib')
    joblib.dump({'model':residual,'feature_names':names,'report':{'target':'ddg','representation':'surface',
        'prior_method':'frozen_prior','method':'prior_residual'}},tmp_path/'residual.joblib')
    bundle_mutation_ensemble(tmp_path/'base.joblib',tmp_path/'residual.joblib',tmp_path/'ensemble.joblib')
    np.savez(tmp_path/'base.npz',X=x[::-1],sample_ids=ids[::-1],feature_names=names)
    np.savez(tmp_path/'residual.npz',X=x,sample_ids=ids,feature_names=names)
    prior = pd.DataFrame({'sample_id':ids,'predicted':np.arange(8)/10.,'target':'ddg',
                          'method':'frozen_prior','status':'ok','error':''})
    expected = .5*(base.predict(x)+residual.predict(prior.predicted,x))
    prior.loc[2,['predicted','status','error']] = [np.nan,'error','unavailable residue']
    expected[2] = base.predict(x[[2]])[0]
    prior = prior.iloc[[4,1,0,7,3,2,6,5]]
    result = predict_mutation_ensemble(tmp_path/'ensemble.joblib',tmp_path/'base.npz',tmp_path/'residual.npz',prior)
    assert result.sample_id.tolist() == prior.sample_id.tolist()
    ordered = result.set_index('sample_id').loc[ids]
    np.testing.assert_allclose(ordered.predicted,expected,atol=1e-12)
    assert ordered.loc['s2','status'] == 'base_only'
    assert ordered.loc['s2','error'] == 'unavailable residue'
    assert ordered.loc['s2','prior_status'] == 'error'
    assert ordered.loc['s0','components_used'] == 'equal_ensemble'


@pytest.mark.parametrize('ids', [['a', 'b', 'c', 'd'], ['001', '002', '003', '004']])
def test_scalar_prior_fusion_preserves_sample_alignment(tmp_path, ids):
    from topoppi_affinity.mutation_ensemble import bundle_scalar_prior_ensemble, predict_scalar_prior_ensemble
    x = np.arange(12, dtype=float).reshape(4, 3)
    ids, names = np.asarray(ids), np.asarray(['plm.a','contact.b','surface.c'])
    model = fit_frozen_regressor(x, np.arange(4.), ids, family='ridge',params={'alpha':1.},antisymmetric=True)
    joblib.dump({'model':model,'feature_names':names,'report':{'target':'ddg'}},tmp_path/'base.joblib')
    bundle_scalar_prior_ensemble(tmp_path/'base.joblib','BA-Cycle-partner-MC8',tmp_path/'ensemble.joblib')
    np.savez(tmp_path/'features.npz',X=x[::-1],sample_ids=ids[::-1],feature_names=names)
    prior = pd.DataFrame({'sample_id':ids,'predicted':[.1,.2,np.nan,.4],'target':'ddg',
                         'method':'BA-Cycle-partner-MC8','status':['ok','ok','error','ok'],'error':['','','missing','']})
    prior.iloc[[3,1,0,2]].to_csv(tmp_path / 'prior.csv', index=False)
    result = predict_scalar_prior_ensemble(tmp_path/'ensemble.joblib',tmp_path/'features.npz',tmp_path/'prior.csv')
    expected = .5*(model.predict(x)+prior.predicted)
    expected[2] = model.predict(x[[2]])[0]
    np.testing.assert_allclose(result.set_index('sample_id').loc[ids,'predicted'],expected)
    assert result.loc[result.sample_id.eq(ids[2]),'components_used'].item() == 'base_only'
