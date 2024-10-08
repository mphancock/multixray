from pathlib import Path

import numpy as np

import sys
sys.path.append("..")
from multi_state_multi_condition_model import MultiStateMultiConditionModel
from charmm import CHARMMRestraint, charmm_restraints


if __name__ == "__main__":
    w_mat = np.array([[1.0]])

    msmc_m = MultiStateMultiConditionModel(
        pdb_files=[pdb_file],
        w_mat=w_mat,
        crystal_symmetries=None
    )

    h, m = msmc_m.get_hs()[0], msmc_m.get_m()

    charmm_rs = charmm_restraints(
        m=m,
        h=h,
        eps=False
    )


