from pathlib import Path
import sys

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import cctbx.miller
import iotbx
from scitbx.array_family import flex

sys.path.append(str(Path(Path.home(), "xray/src")))
import miller_ops


def get_f_model(
        pdb_file,
        uc_dim,
        sg_symbol,
        res
):
    crystal_symmetry = cctbx.xray.crystal.symmetry(
        unit_cell=uc_dim,
        space_group_symbol=sg_symbol
    )

    xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    # miller_set = cctbx.miller.build_set(
    #   crystal_symmetry=crystal_symmetry,
    #   anomalous_flag=False,
    #   d_min=2.0,
    # )

    f_model = xray_structure.structure_factors(
        d_min=1.0
    ).f_calc().amplitudes()
    print(f_model.size())

    return f_model


def get_f_model_from_cif(
        pdb_file,
        cif_file,
        res=None
):
    f_obs = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.F_meas_au"
    )
    # f_obs = miller_ops.filter_f_obs_resolution(
    #     f_obs=f_obs,
    #     d_min=res,
    #     d_max=None
    # )
    crystal_symmetry = f_obs.crystal_symmetry()

    if res:
        f_obs = f_obs.complete_array(
            d_min=res,
            d_max=None
        )

        f_obs.generate_r_free_flags(
            fraction=.1
        )


    xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        target_name="ls"
    )

    f_model = f_model_manager.f_model()

    print(f_model.size())

    return f_model


if __name__ == "__main__":
    pdb_dir = Path(Path.home(), "xray/dev/19_synthetic_native_2/data/pdbs/4_state_2")
    pdb_files = list(pdb_dir.glob("*.pdb"))

    for pdb_file in pdb_files:
        print(pdb_file)
        model_cif_file = Path(Path.home(), "xray/dev/19_synthetic_native_2/data/cifs/{}/{}.cif".format(pdb_dir.name, pdb_file.stem))

        f_obs_array = get_f_model(
            pdb_file=pdb_file,
            uc_dim=(58.3050, 36.1540, 25.3620, 90.0000, 103.0900, 90.0000),
            sg_symbol="C 1 2 1",
            res=1
        )

        flags_array = f_obs_array.generate_r_free_flags(
            fraction=.1
        )

        status = flex.std_string()

        for i in range(flags_array.size()):
            if flags_array.data()[i]:
                status.append("f")
            else:
                status.append("o")

        status_array = flags_array.customized_copy(data=status)
        print(status_array.data())
        print(flags_array.data())

        # for i in range(flags_array.size()):
        #     if flags_array.data()[i]:
        #         status_array.data()[i] = "f"
        #     else:
        #         status_array.data()[i] = "o"


        # # f_model = get_f_model_from_cif(
        # #     pdb_file=pdb_file,
        # #     cif_file=cif_file,
        # #     res=1.0
        # # )

        # # Set flags.
        # # status_array = miller_ops.get_miller_array(
        # #     f_obs_file=cif_file,
        # #     label="_refln.status"
        # # )
        # # flags_array = status_array.customized_copy(data=status_array.data()=="f")
        # # f_obs_array, flags_array = f_model.common_sets(other=flags_array)

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
