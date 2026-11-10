import numpy as np
from sklearn.ensemble import ExtraTreesRegressor

from topoppi_affinity.legacy_antisymmetry import ConditionalOddTree, reverse_legacy


def test_reverse_involution_fixed_context_and_odd_projection():
    names = np.array(['sequence.aa_A','mutation.wt_A','mutation.mt_A','mutation.delta_A',
                      'mutation.contact_delta_charge_degree','mutation.absdelta_charge',
                      'epistasis.delta_charge_product_proximity'])
    x = np.random.default_rng(11).normal(size=(12,len(names)))
    reverse = reverse_legacy(x,names)
    np.testing.assert_array_equal(reverse_legacy(reverse,names),x)
    np.testing.assert_array_equal(reverse[:,[0,5,6]],x[:,[0,5,6]])
    model = ConditionalOddTree(ExtraTreesRegressor(n_estimators=4,random_state=2).fit(x,np.arange(len(x))),names)
    np.testing.assert_allclose(model.predict(x),-model.predict(reverse),atol=0,rtol=0)
    noop = x.copy()
    noop[:,2]=noop[:,1]
    noop[:,[3,4]]=0
    np.testing.assert_array_equal(model.predict(noop),np.zeros(len(noop)))
