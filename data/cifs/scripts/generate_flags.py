from pathlib import Path
import random
random.seed(0)
import sys

import cctbx
import cctbx.miller
import iotbx
import iotbx.cif
import iotbx.reflection_file_reader as cif_input

sys.path.append(str(Path(Path.home(), "xray/src")))
import miller_ops


if __name__ == "__main__":
    f_obs_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")
    status_array = miller_ops.get_miller_array(
        f_obs_file=f_obs_file,
        label="_refln.status"
    )

    f_obs_array = miller_ops.get_miller_array(
        f_obs_file=f_obs_file,
        label="_refln.F_meas_au"
    )
    flags_array = status_array.customized_copy(data=status_array.data()=="f")

    print(f_obs_array.size(), flags_array.size())
    f_obs_array = miller_ops.clean_miller_array(f_obs_array)
    flags_array = miller_ops.clean_miller_array(flags_array)
    print(f_obs_array.size(), flags_array.size())

    f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

    for i in range(f_obs_array.size()):
        print(i, f_obs_array.indices()[i], flags_array.data()[i], f_obs_array.data()[i])

        if i > 50:
            break

    f_work = f_obs_array.select(~(flags_array.data()))
    print(f_work.size())