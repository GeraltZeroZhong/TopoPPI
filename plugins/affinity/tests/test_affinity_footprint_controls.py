import numpy as np
import torch
from torch_geometric.data import Data

from topoppi_affinity.footprint_controls import permute_footprints


def test_only_residue_assignment_changes_and_both_states_share_permutation():
    x = torch.arange(12*45,dtype=torch.float32).reshape(12,45)
    data = Data(x=x,x_mut=x.clone(),side=torch.tensor([0]*6+[1]*6),
                contact_mask=torch.tensor([True]*3+[False]*3+[True]*3+[False]*3))
    result = permute_footprints(data,11)
    torch.testing.assert_close(result.x[:,:30],data.x[:,:30],atol=0,rtol=0)
    torch.testing.assert_close(result.x,result.x_mut,atol=0,rtol=0)
    assert not torch.equal(result.x,data.x)
    for side in [0,1]:
        for contact in [False,True]:
            mask = (data.side==side)&(data.contact_mask==contact)
            np.testing.assert_array_equal(np.sort(result.x[mask,30:].numpy(),axis=0),
                                          np.sort(data.x[mask,30:].numpy(),axis=0))
