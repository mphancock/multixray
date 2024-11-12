from pathlib import Path

import numpy as np

import sys
sys.path.append("..")
from multi_state_multi_condition_model import MultiStateMultiConditionModel
from align_imp import get_multi_state_multi_cond_rmsd, compute_rmsd_between_ordered_states


if __name__ == "__main__":
    msmc_m_0 = MultiStateMultiConditionModel(
        pdb_files=[Path("./pdbs/test_0.pdb")],
        w_mat=np.array([[0.5,0.5],[0.5,0.5]]),
        crystal_symmetries=None
    )

    msmc_m_1 = MultiStateMultiConditionModel(
        pdb_files=[Path("./pdbs/test_1.pdb")],
        w_mat=np.array([[0.5,0.75],[0.5,0.25]]),
        crystal_symmetries=None
    )
    ## avg positions of model 1 is
    ## 0.5, 0.5, 0.5
    ## 0.75*1+0.25*0, 0.75*0+0.25*1, 0.75*0+0.25*1 = 0.75, 0.25, 0.25

    ## sqrt ( 1/1 * (0.5**2 + 0.5**2 + 0.5**2)) = 0.86602540378
    print(get_multi_state_multi_cond_rmsd(msmc_m_0, msmc_m_1, cond=0))

    ## sqrt ( 1/1 * ((0-0.75)**2 + (0-0.25)**2 + (0-0.25)**2)) = 0.82915619758
    print(get_multi_state_multi_cond_rmsd(msmc_m_0, msmc_m_1, cond=1))

    print(get_multi_state_multi_cond_rmsd(msmc_m_0, msmc_m_1, cond=None))

    ## TEST WITH ACUTAL PDB STRUCTURES
    msmc_m_2 = MultiStateMultiConditionModel(
        pdb_files=[Path("./pdbs/test_2.pdb")],
        w_mat=np.array([[1.0]]),
        crystal_symmetries=None
    )

    msmc_m_3 = MultiStateMultiConditionModel(
        pdb_files=[Path("./pdbs/test_3.pdb")],
        w_mat=np.array([[1.0]]),
        crystal_symmetries=None
    )

    ## 0.297
    print(get_multi_state_multi_cond_rmsd(msmc_m_2, msmc_m_3, cond=None))

    ## TEST ORDERED
    ## model 0 state 0 is 0 0 0
    ## model 0 state 1 is 0 0 0
    ## model 1 state 0 is 1 0 0
    ## model 1 state 1 is 0 1 1
    ## rmsd between model 0 state 0 and model 1 state 0 is
    ## sqrt( 1/1 * (1**2 + 0**2 + 0**2)) = 1
    ## rmsd between model 0 state 1 and model 1 state 1 is
    ## sqrt( 1/1 * (1**2 + 1**2 + 1**2)) = 1.41421356237
    ## avg rmsd is 1.20710678118
    print(compute_rmsd_between_ordered_states(msmc_m_0, msmc_m_1, cond=0))