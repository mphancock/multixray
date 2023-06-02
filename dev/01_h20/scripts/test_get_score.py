from pathlib import Path
import sys

import IMP
import IMP.atom

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import iotbx

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
import miller_ops
import xray_struct
import cctbx_score


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean_h20.pdb")
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    cif_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")
    f_obs_array = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.F_meas_au"
    )
    f_obs_array = miller_ops.clean_miller_array(f_obs_array)

    # Set flags.
    status_array = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.status"
    )
    flags_array = status_array.customized_copy(data=status_array.data()=="f")
    f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

    score_dict = cctbx_score.get_score(
        m=m,
        uc_dim=(58.305,36.154,25.362,90.00,103.09,90.00),
        sg_symbol="C 1 2 1",
        f_obs=f_obs_array,
        r_free_flags=flags_array,
        target="ml"
    )
    # print(score_dict)