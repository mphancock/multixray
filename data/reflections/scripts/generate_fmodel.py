from pathlib import Path
import sys

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import iotbx

sys.path.append(str(Path(Path.home(), "xray/src")))
import miller_ops


def get_f_model(
        pdb_file,
        cif_file,
        res
):
    f_obs = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.F_meas_au"
    )
    f_obs = miller_ops.filter_f_obs_resolution(
        f_obs=f_obs,
        d_min=res,
        d_max=None
    )
    crystal_symmetry = f_obs.crystal_symmetry()

    xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        target_name="ls"
    )

    f_model = f_model_manager.f_model()

    return f_model


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhf_refine.pdb")

    cif_file = Path(Path.home(), "xray/data/reflections/7mhf/7mhf.cif")
    model_cif_file = Path(Path.home(), "xray/data/reflections/7mhf/7mhf_refine.cif")

    f_model = get_f_model(
        pdb_file=pdb_file,
        cif_file=cif_file,
        res=None
    )

    # Set flags.
    status_array = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.status"
    )
    flags_array = status_array.customized_copy(data=status_array.data()=="f")
    f_obs_array, flags_array = f_model.common_sets(other=flags_array)

    cif_model = iotbx.cif.model.cif()

    cif_block = iotbx.cif.miller_arrays_as_cif_block(
        array=f_obs_array,
        # array_type="meas",
        column_names=['_refln.F_meas_au','_refln.F_meas_sigma_au'],
        miller_index_prefix="_refln",
        format="mmcif"
    )
    cif_block.add_miller_array(
        array=status_array,
        column_name="_refln.status"
    )

    # print(type(cif_block.cif_block))
    cif_model["Global"] = cif_block.cif_block

    with open(str(model_cif_file), "w") as f:
      print(cif_model, file=f)
