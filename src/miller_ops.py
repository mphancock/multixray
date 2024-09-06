from pathlib import Path

import iotbx.reflection_file_reader as cif_input


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


def clean_miller_array(
        miller_array
):
    miller_array = miller_array.map_to_asu()
    miller_array = miller_array.merge_equivalents().array()
    miller_array = miller_array.remove_systematic_absences()

    return miller_array


def filter_f_obs_resolution(
        f_obs,
        d_max,
        d_min
):
    sel = f_obs.resolution_filter_selection(
        d_max=d_max,
        d_min=d_min
    )
    f_obs_filter = f_obs.select(
        selection=sel
    )

    return f_obs_filter


