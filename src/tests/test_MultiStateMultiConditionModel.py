from pathlib import Path
import numpy as np

import IMP
import IMP.atom

import sys
sys.path.append("..")
from multi_state_multi_condition_model import MultiStateMultiConditionModel


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "Documents/xray/data/pdbs/3k0m/3k0m_refine.pdb")
    w_mat = np.array([[0.5, 1.0],[0.5, 0.0]])
    msmc_model = MultiStateMultiConditionModel(pdb_files=[pdb_file, pdb_file], w_mat=w_mat)

    ## check all occupancies
    multi_xray_struct = msmc_model.get_multi_xray_structure(0)
    print(multi_xray_struct.scatterers()[0].occupancy)
    print(multi_xray_struct.scatterers()[-1].occupancy)

    multi_xray_struct = msmc_model.get_multi_xray_structure(1)
    print(multi_xray_struct.scatterers()[0].occupancy)
    print(multi_xray_struct.scatterers()[-1].occupancy)

    ## change occupancies
    msmc_model.set_occs_for_condition_i([0.25, 0.75], 0)
    multi_xray_struct = msmc_model.get_multi_xray_structure(0)
    print(multi_xray_struct.scatterers()[0].occupancy)
    print(multi_xray_struct.scatterers()[-1].occupancy)

    multi_xray_struct = msmc_model.get_multi_xray_structure(0)
    print(multi_xray_struct.scatterers()[0].site)

    pid = msmc_model.get_pids()[0]
    d = IMP.core.XYZR(msmc_model.get_m(), pid)
    d.set_coordinates(IMP.algebra.Vector3D(0, 0, 0))

    multi_xray_struct = msmc_model.get_multi_xray_structure(0)
    print(multi_xray_struct.scatterers()[0].site)

    pdb_0_file = Path(Path.home(), "Documents/xray/tmp/test_0.pdb")
    pdb_1_file = Path(Path.home(), "Documents/xray/tmp/test_1.pdb")
    msmc_model.write_state_pdb_file(0, pdb_0_file)
    msmc_model.write_state_pdb_file(1, pdb_1_file)

    multi_pdb_file = Path(Path.home(), "Documents/xray/tmp/test.pdb")
    msmc_model.write_pdb_file(multi_pdb_file)

    msmc_model_2 = MultiStateMultiConditionModel(pdb_files=[pdb_0_file, pdb_1_file], w_mat=w_mat)
    msmc_model_3 = MultiStateMultiConditionModel(pdb_files=[multi_pdb_file], w_mat=w_mat)