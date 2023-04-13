from pathlib import Path
import sys

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import iotbx

sys.path.append("/home/matthew/xtal_python/src")
import xray_struct
import miller_ops


"""
Deprecated. 
"""
def get_r_factor(
        pdb_file,
        uc_dim,
        sg_symbol,
        cif_file
):
    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=uc_dim,
        space_group_symbol=sg_symbol
    )

    xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    f_obs = miller_ops.get_f_obs(str(cif_file))

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        target_name="ml"
    )

    f_model_manager.update_all_scales(
        remove_outliers=False
    )

    r_factor = f_model_manager.r_all()

    return r_factor


"""
Given a set of IMP hierarchies, the unit cell dimensions, the space group symbol, the cif file, and the target function name; the function returns the associated score. The target functions include "rf", "ls", and "ml".

flag_type: specifies what set of miller indices to evaluate the scores over. 
"""
# def get_score_from_cif_file(
#         m,
#         uc_dim,
#         sg_symbol,
#         cif_file,
#         res,
#         target_name,
#         flag_type=None
# ):
#
#
#
#     get_score(
#         m=m,
#         uc_dim=uc_dim,
#         sg_symbol=sg_symbol,
#         f_obs_work
#     )


def get_score(
        m,
        uc_dim,
        sg_symbol,
        f_obs_work,
        target_name
):
    # print(f_obs_work.size())

    r_factor_target = False
    if target_name == "rf":
        r_factor_target = True

    # xray_structure = get_cctbx_multi_structure(
    #     pdb_files=pdb_files,
    #     weights=weights,
    #     crystal_symmetry=crystal_symmetry
    # )

    xray_structure = xray_struct.get_xray_structure(
        m=m,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol
    )

    xray_structure.scatterers().flags_set_grads(
        state=False
    )
    xray_structure.scatterers().flags_set_grad_site(
        iselection=xray_structure.all_selection().iselection()
    )

    if r_factor_target == "rf":
        target = "ml"
    else:
        target = target_name

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs_work,
        target_name=target
    )
    f_model_manager.update_all_scales()

    if r_factor_target:
        r_factor = f_model_manager.r_all()
        #
        # print("r_all: {}".format(f_model_manager.r_all()))
        # print("r_work: {}".format(f_model_manager.r_work()))
        # print("r_free: {}".format(f_model_manager.r_free()))

        return r_factor, None

    fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)
    fmodels.update_xray_structure(
        xray_structure=xray_structure,
        update_f_calc=True)

    fmodels_target_and_gradients = fmodels.target_and_gradients(
        compute_gradients=True)
    score = fmodels_target_and_gradients.target()
    grads = fmodels_target_and_gradients.gradients()

    return score, grads


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    pdb_files = [pdb_file]

    cif_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")

    uc_dim = "58.3050, 36.1540, 25.3620, 90.0000, 103.0900, 90.0000"
    sg_symbol = "C 1 2 1"

    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=uc_dim,
        space_group_symbol=sg_symbol
    )

    # xray_structure.show_scatterers()
    # crystal_symmetry = cctbx.crystal.symmetry(
    #     unit_cell=uc_dim,
    #     space_group_symbol=sg_symbol
    # )

    xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    f_obs = miller_ops.get_f_obs(str(cif_file))
    f_obs.show_comprehensive_summary()
    print(f_obs.data()[0])

    f_c = f_obs.copy()
    f_c.data()[0] = 0
    print(f_c.data()[0])
    f_c = f_c.structure_factors_from_scatterers(
        algorithm="fft",
        xray_structure=xray_structure
    ).f_calc()
    print(f_c.as_amplitude_array().data()[0])

    # dir = cctbx.xray.structure_factors.from_scatterers_direct(
    #     xray_structure,
    #     f_obs
    # )
    # dir_f_calc = dir.f_calc()
    # # print(type(dir_f_calc))
    #
    # print(dir_f_calc.as_amplitude_array().data()[0])
    xray_structure = xray_struct.get_cctbx_multi_structure_from_pdbs(
        pdb_files=pdb_files,
        weights=[1.],
        crystal_symmetry=crystal_symmetry
    )

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        target_name="ls"
    )

    f_model_manager.update_all_scales(
        remove_outliers=False
    )

    print(f_model_manager.f_calc().as_amplitude_array().data()[0])

    f_model_manager.show()
    print(f_model_manager.r_all())

