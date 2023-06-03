from pathlib import Path
import sys

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import iotbx

sys.path.append("/home/matthew/xtal_python/src")
import xray_struct
import miller_ops


def get_score(
        m,
        f_obs,
        r_free_flags,
        target
):
    crystal_symmetry = f_obs.crystal_symmetry()
    xray_structure = xray_struct.get_xray_structure(
        m=m,
        crystal_symmetry=crystal_symmetry
    )

    xray_structure.scatterers().flags_set_grads(
        state=False
    )
    xray_structure.scatterers().flags_set_grad_site(
        iselection=xray_structure.all_selection().iselection()
    )
    xray_structure.scatterers().flags_set_grad_occupancy(
        iselection=xray_structure.all_selection().iselection()
    )

    r_factor_target = False
    target_name = target

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        target_name=target_name
    )
    f_model_manager.update_all_scales()

    r_work = f_model_manager.r_work()
    r_free = f_model_manager.r_free()
    r_all = f_model_manager.r_all()

    fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)
    fmodels.update_xray_structure(
        xray_structure=xray_structure,
        update_f_calc=True)

    fmodels_target_and_gradients = fmodels.target_and_gradients(
        compute_gradients=True)
    score = fmodels_target_and_gradients.target()
    grads = fmodels_target_and_gradients.gradients()

    site_grads = list()
    occ_grads = list()

    for i in range(len(xray_structure.scatterers())):
        site_grads.append(grads[i*4:i*4+3])
        occ_grads.append(grads[3+i*4])

    results_dict = dict()
    results_dict["score"] = score
    results_dict["grads_site"] = site_grads
    results_dict["grads_occ"] = occ_grads
    results_dict["r_work"] = r_work
    results_dict["r_free"] = r_free
    results_dict["r_all"] = r_all

    return results_dict