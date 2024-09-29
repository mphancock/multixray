from pathlib import Path
import numpy as np

import IMP
import IMP.atom

import sys
sys.path.append("..")
from multi_state_multi_condition_model import MultiStateMultiConditionModel
from xray_restraint import XtalRestraint
from miller_ops import get_miller_array, clean_miller_array


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "Documents/xray/data/pdbs/3k0m/3k0m.pdb")
    w_mat = np.array([[0.5, 1.0],[0.5, 0.0]])
    msmc_m = MultiStateMultiConditionModel(
        pdb_file=pdb_file,
        w_mat=w_mat
    )

    cif_file = Path(Path.home(), "Documents/xray/data/cifs/3k0m/3k0m.cif")
    f_obs_array = get_miller_array(
        f_obs_file=cif_file,
        # label="_refln.F_meas_au"
        label="_refln.intensity_meas"
    )
    f_obs_array = clean_miller_array(f_obs_array)

    # Set flags from file.
    status_array = get_miller_array(
        f_obs_file=cif_file,
        label="_refln.status"
    )
    flags_array = status_array.customized_copy(data=status_array.data()=="f")
    f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

    pid = msmc_m.get_all_pids()[0]
    d = IMP.core.XYZR(msmc_m.get_m(), pid)

    print(d.get_derivatives())

    r_xray = XtalRestraint(
        msmc_m=msmc_m,
        cond=0,
        f_obs=f_obs_array,
        free_flags=flags_array,
        w_xray=1,
        update_scale=True,
        update_k1=True,
        update_freq=1,
        ref_com=None
    )
    r_xray.evaluate(True)

    print(d.get_derivatives())