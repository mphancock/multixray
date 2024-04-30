from pathlib import Path
import sys
import numpy as np

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
        res,
        ws=None
):
    crystal_symmetry = cctbx.xray.crystal.symmetry(
        unit_cell=uc_dim,
        space_group_symbol=sg_symbol
    )

    xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    if ws:
        n_scatt = xray_structure.scatterers().size()
        n_scatt_per_state = int(n_scatt/len(ws))

        for i in range(len(ws)):
            for j in range(n_scatt_per_state):
                xray_structure.scatterers()[i*n_scatt_per_state+j].occupancy = ws[i]

    # xray_structure.show_scatterers()

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


def get_f_model_from_f_obs(
        pdb_file,
        cif_file,
        occs,
        res=None
):
    f_obs = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.F_meas_au"
    )

    f_obs = miller_ops.clean_miller_array(f_obs)
    status_array = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.status"
    )
    flags = status_array.customized_copy(data=status_array.data()=="f")
    f_obs, flags = f_obs.common_sets(other=flags)

    f_obs = miller_ops.filter_f_obs_resolution(
        f_obs=f_obs,
        d_max=None,
        d_min=None
    )
    flags = miller_ops.filter_f_obs_resolution(
        f_obs=flags,
        d_max=None,
        d_min=None
    )
    crystal_symmetry = f_obs.crystal_symmetry()

    # if res:
    #     f_obs = f_obs.complete_array(
    #         d_min=res,
    #         d_max=None
    #     )

    #     f_obs.generate_r_free_flags(
    #         fraction=.1
    #     )

    xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    n_scatt = int(xray_structure.scatterers().size())
    n_state = len(occs)
    n_scatt_per_state = n_scatt//n_state

    for i in range(n_scatt_per_state):
        for state in range(n_state):
            xray_structure.scatterers()[i+(state*n_scatt_per_state)].occupancy = occs[state]

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        r_free_flags=flags,
        target_name="ml"
    )
    f_model_manager.update_all_scales(apply_scale_k1_to_f_obs=False,remove_outliers=False)

    f_model = f_model_manager.f_model().as_amplitude_array()

    print(f_model.size())
    print(f_model_manager.r_free())

    return f_model, status_array


def get_status_array(
        flags_array
):
    status = flex.std_string()

    for i in range(flags_array.size()):
        if flags_array.data()[i]:
            status.append("f")
        else:
            status.append("o")

    status_array = flags_array.customized_copy(data=status)

    return status_array


def write_cif(
    f_obs,
    status_array,
    cif_file
):
    cif_model = iotbx.cif.model.cif()

    cif_block = iotbx.cif.miller_arrays_as_cif_block(
        array=f_obs,
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

    with open(str(cif_file), "w") as f:
        print(cif_model, file=f)


def randomize_amplitude(
        f_obs
):
    for i in range(f_obs.size()):
        amp = f_obs.data()[i]

        print(amp)

        mean = amp
        std = amp*.1
        f_obs.data()[i] = np.random.normal(
            loc=mean,
            scale=std
        )

    return f_obs


if __name__ == "__main__":
    pdb_file = Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30/0.pdb")
    cif_file = Path("/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhf.cif")
    out_cif_file = Path(Path.home(), "xray/tmp/1.cif")

    f_model, status_array = get_f_model_from_f_obs(
        pdb_file=pdb_file,
        cif_file=cif_file,
        occs=[0.75,0.25],
    )

    # f_obs_array = randomize_amplitude(
    #     f_obs=f_model
    # )

    # # flags_array = f_obs_array.generate_r_free_flags(
    # #     fraction=.1,
    # #     max_free=None
    # # )

    # for i in range(50):
    #     print(f_obs_array.indices()[i], f_obs_array.data()[i], flags_array.data()[i])

    # status_array = get_status_array(
    #     flags_array=flags_array
    # )

    # print(status_array.data())
    # print(flags_array.data())

    write_cif(
        f_obs=f_model,
        status_array=status_array,
        cif_file=out_cif_file
    )

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


