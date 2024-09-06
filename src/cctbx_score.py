from pathlib import Path
import sys
import time
import pickle

import IMP
import IMP.atom

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray

sys.path.append("/home/matthew/xtal_python/src")
import xray_struct
import miller_ops


def get_score(
        hs,
        occs,
        pids,
        f_obs,
        r_free_flags,
        target,
        ab_file=None,
        update_scale=True,
        update_k1=False,
        u_aniso_file=None,
        delta=None
):
    # print(occs)
    # pid = IMP.atom.Selection(hs[0]).get_selected_particle_indexes()[0]
    # print(IMP.atom.Atom(hs[0].get_model(), pid).get_occupancy())

    # print(f_obs.size())
    crystal_symmetry = f_obs.crystal_symmetry()
    xray_structure = xray_struct.get_xray_structure(
        hs=hs,
        occs=occs,
        pids=pids,
        crystal_symmetry=crystal_symmetry,
        u_aniso_file=u_aniso_file,
        delta=delta
    )

    # print(xray_structure.scatterers()[0].occupancy)
    # print(xray_structure.show_summary())

    xray_structure.scatterers().flags_set_grads(
        state=False
    )
    xray_structure.scatterers().flags_set_grad_site(
        iselection=xray_structure.all_selection().iselection()
    )
    xray_structure.scatterers().flags_set_grad_occupancy(
        iselection=xray_structure.all_selection().iselection()
    )

    target_name = target

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        target_name=target_name
    )

    if update_scale:
        f_model_manager.update_all_scales(
            apply_scale_k1_to_f_obs=update_k1,
            remove_outliers=False
        )

    if ab_file:
        f = open(ab_file, 'rb')
        alpha_beta = pickle.load(f)
        f_model_manager.alpha_beta_cache = alpha_beta

    # f_model_manager.show()

    r_work = f_model_manager.r_work()
    r_free = f_model_manager.r_free()
    r_all = f_model_manager.r_all()

    # xray_structure.show_scatterers()
    # print(len(xray_structure.scatterers()))
    # scatt = xray_structure.scatterers()[0]
    # scatt.show()

    fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)
    fmodels.update_xray_structure(
        xray_structure=xray_structure,
        update_f_calc=True
    )

    fmodels_target_and_gradients = fmodels.target_and_gradients(compute_gradients=True)
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