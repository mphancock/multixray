from pathlib import Path

import iotbx.reflection_file_reader as cif_input


if __name__ == "__main__":
    cif_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")

    reader = cif_input.any_reflection_file(file_name=str(cif_file))
    miller_arrays = reader.as_miller_arrays()

    for miller_array in miller_arrays:
        print(miller_array.crystal_symmetry())
        break

    # flags = None
    # for miller_array in miller_arrays:
    #     miller_array_filt = miller_array.remove_systematic_absences()
    #     if miller_array_filt.is_string_array():
    #         flags = miller_array_filt

    # for i in range(flags.size()):
    #     print(flags.data()[i])

