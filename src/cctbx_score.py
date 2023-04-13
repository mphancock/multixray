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
        uc_dim,
        sg_symbol,
        f_obs,
        r_free_flags,
        target
):
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

    r_factor_target = False
    if target in ["r_work", "r_all", "r_free"]:
        target_name = "ml"
        r_factor_target = True
    else:
        target_name = target

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        target_name=target_name
    )
    f_model_manager.update_all_scales()

    if r_factor_target:
        if target == "r_work":
            r_factor = f_model_manager.r_work()
        elif target == "r_free":
            r_factor = f_model_manager.r_free()
        else:
            r_factor = f_model_manager.r_all()

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
