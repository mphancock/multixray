
from pathlib import Path
import sys

import iotbx.reflection_file_reader as cif_input
from mmtbx.tls import tools
import mmtbx.model
import iotbx.pdb
from cctbx.array_family import flex
from mmtbx.solvent import ensemble_ordered_solvent


def get_miller_array(
        f_obs_file,
        label
):
    reader = cif_input.any_reflection_file(file_name=str(f_obs_file))
    miller_arrays = reader.as_miller_arrays()

    array_id = -1
    for i in range(len(miller_arrays)):
        miller_array = miller_arrays[i]
        # if '_refln.F_meas_au' in miller_array.info().labels:
        if label in miller_array.info().labels:
            array_id = i

    label_array = miller_arrays[array_id]

    if array_id < 0:
        raise RuntimeError("{} was not found in {}".format(label, f_obs_file))

    return label_array


if __name__ == "__main__":
  pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhj_heavy.pdb")
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(str(pdb_file)))
  xray_struct = model.get_xray_structure()

  f_obs_file = Path(Path.home(), "xray/data/cifs/7mhf/7mhj.cif")
  f_obs_array = get_miller_array(
      f_obs_file=f_obs_file,
      label="_refln.F_meas_au"
  )
  status_array =get_miller_array(
      f_obs_file=f_obs_file,
      label="_refln.status"
  )
  flags_array = status_array.customized_copy(data=status_array.data()=="f")
  f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg_params.algorithm = "direct"
  sfg_params.cos_sin_table = False

  # params = iotbx.phil.parse(master_params_str).extract()
  # print(params.preserved_solvent_minimum_distance)

  fmodel = mmtbx.f_model.manager(
    xray_structure    = xray_struct,
    f_obs             = f_obs_array,
    r_free_flags      = flags_array,
    target_name       = "ls_wunit_k1",
    sf_and_grads_accuracy_params = sfg_params
  )

  manager = ensemble_ordered_solvent.manager(
      fmodel=fmodel,
      model=model,
      verbose=1
  )