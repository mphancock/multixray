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
sys.path.append(str(Path(Path.home(), "Documents/xray/src")))
import miller_ops


def get_f_model(
        pdb_file,
        uc_dim,
        sg_symbol,
        res,
        ws=None,

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

    # f_model = xray_structure.structure_factors(
    #     d_min=res
    # ).f_calc().amplitudes()
    # print(f_model.d_min())

    f_c = abs(xray_structure.structure_factors(d_min=res).f_calc())
    r_free_flags = f_c.generate_r_free_flags(fraction = 0.1, max_free = 99999999)

    mask_params = mmtbx.masks.mask_master_params.extract()
    fmodel = mmtbx.f_model.manager(
        # mask_params    = mask_params,
        r_free_flags   = r_free_flags,
        f_obs          = f_c,
        xray_structure = xray_structure)
    fmodel.update_xray_structure(
        xray_structure = xray_structure,
        update_f_calc  = True,
        update_f_mask  = True)
    fmodel_kbu = fmodel.fmodel_kbu()
    fmodel_kbu.update(
        k_sols = 1,
        b_sols = 50,
        b_cart = [0,0,0,0,0,0])

    return abs(fmodel_kbu.f_model)


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
        f_obs,
        dist="norm",
        std=None
):
    for i in range(f_obs.size()):
        amp = f_obs.data()[i]

        # print(amp)

        mean = amp

        if dist == "norm":
            f_obs.data()[i] = np.random.normal(
                loc=mean,
                scale=std
            )
        elif dist == "uni":
            if mean-std < 0:
                low = 0
            else:
                low = mean-std

            f_obs.data()[i] = np.random.uniform(
                low=low,
                high=mean+std
            )


    return f_obs


def get_f_obs_freer(
    d_min,
    k_sol,
    b_sol,
    b_cart,
    xray_structure,
    radial_shell_width=None
):
  f_dummy = abs(xray_structure.structure_factors(d_min = d_min,
    anomalous_flag = False).f_calc())
  r_free_flags = f_dummy.generate_r_free_flags(fraction = 0.1,
                                               max_free = 99999999)
  mask_params = mmtbx.masks.mask_master_params.extract()
  if( radial_shell_width is not None ):
    mask_params.radial_shell_width = radial_shell_width
  if( type(k_sol) is list):
    mask_params.n_radial_shells = len(k_sol)
  fmodel = mmtbx.f_model.manager(
    mask_params    = mask_params,
    r_free_flags   = r_free_flags,
    f_obs          = f_dummy,
    xray_structure = xray_structure)
  fmodel.update_xray_structure(
    xray_structure = xray_structure,
    update_f_calc  = True,
    update_f_mask  = True)
  fmodel_kbu = fmodel.fmodel_kbu()
  fmodel_kbu.update(
    k_sols = k_sol,
    b_sols = b_sol,
    b_cart = b_cart)
  f_obs = abs(fmodel_kbu.f_model)
  return f_obs, r_free_flags


