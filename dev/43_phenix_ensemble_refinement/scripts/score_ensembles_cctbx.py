from pathlib import Path
import pickle

import IMP
import IMP.atom

import mmtbx.f_model
import mmtbx.model
import cctbx.crystal
import cctbx.xray
import iotbx

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
sys.path.append(str(Path(Path.home(), "xray/data/cifs/scripts")))
# sys.path.append("../src")
from miller_ops import get_miller_array, clean_miller_array


if __name__ == "__main__":
    score_dir = Path(Path.home(), "xray/dev/32_score")

    pdb_file = Path("../data/pdbs/3k0m.pdb")
    f_obs_file = Path(Path.home(), "xray/data/cifs/3k0m/3k0m.cif")

    f_obs = get_miller_array(
            f_obs_file=f_obs_file,
            label="_refln.intensity_meas"
    )
    print(f_obs.indices()[0], f_obs.data()[0])

    f_obs = clean_miller_array(f_obs)
    status_array = get_miller_array(
        f_obs_file=f_obs_file,
        label="_refln.status"
    )
    flags = status_array.customized_copy(data=status_array.data()=="f")
    f_obs, flags = f_obs.common_sets(other=flags)

    # d_min = 2.0
    # f_obs = miller_ops.filter_f_obs_resolution(
    #         f_obs=f_obs,
    #         d_max=None,
    #         d_min=d_min
    # )
    # flags = miller_ops.filter_f_obs_resolution(
    #         f_obs=flags,
    #         d_max=None,
    #         d_min=d_min
    # )

    crystal_symmetry = f_obs.crystal_symmetry()

    model = mmtbx.model.manager(
        model_input=iotbx.pdb.input(str(pdb_file)),
        crystal_symmetry=crystal_symmetry
    )

    xray_structure = model.get_xray_structure()
    print(xray_structure.scatterers().size())

    # m = IMP.Model()
    # hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    # pids = list()
    # for h in hs:
    #     pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

    # score_result = cctbx_score.get_score(
    #     hs=hs,
    #     occs=occs,
    #     pids=pids,
    #     f_obs=f_obs,
    #     r_free_flags=flags,
    #     target="ml",
    #     update_k1=True
    # )
    # print(score_result["score"], score_result["r_work"], score_result["r_free"])

    # xray_structure = xray_struct.get_xray_structure(
    #     hs=hs,
    #     occs=occs,
    #     pids=pids,
    #     crystal_symmetry=crystal_symmetry
    # )

    # n_scatt = int(xray_structure.scatterers().size())
    # print(type(n_scatt))
    # n_scatt_state = n_scatt//len(occs)
    # n_state = len(occs)

    # for i in range(n_scatt_state):
    #     for j in range(n_state):
    #         xray_structure.scatterers()[n_scatt_state*j+i].occupancy = occs[j]

    # print(xray_structure.scatterers()[0].occupancy)
    # print(xray_structure.scatterers()[n_scatt_state-1].occupancy)
    # print(xray_structure.scatterers()[-1].occupancy)

    # xray_structure.scatterers().flags_set_grads(
    #     state=False
    # )
    # xray_structure.scatterers().flags_set_grad_site(
    #     iselection=xray_structure.all_selection().iselection()
    # )
    # xray_structure.scatterers().flags_set_grad_occupancy(
    #     iselection=xray_structure.all_selection().iselection()
    # )

    # for i in range(18704):
    #     print(i, flags.data()[i], f_obs.data()[i])

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        r_free_flags=flags,
        target_name="ml",
        max_number_of_bins=1
    )
    f_model_manager.update_all_scales(apply_scale_k1_to_f_obs=True,remove_outliers=True)

    r_work = f_model_manager.r_work()
    r_free = f_model_manager.r_free()
    r_all = f_model_manager.r_all()

    fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)
    # fmodels.update_xray_structure(
    #     xray_structure=xray_structure,
    #     update_f_calc=True
    # )

    fmodels_target_and_gradients = fmodels.target_and_gradients(compute_gradients=True)
    score = fmodels_target_and_gradients.target()
    grads = fmodels_target_and_gradients.gradients()

    # xray_structure.show_scatterers()
    # scatt = xray_structure.scatterers()[0]
    # # print(dir(scatt))
    # print(len(xray_structure.scatterers()))
    # scatt.show()
    print(r_free, flags.size(), r_work, f_obs.size())
    print("likelihood:", score)
    print("r_free:", r_free)
    print("r_work:", r_work)

    # for i in range(len(occs)):
    #     for pid in IMP.atom.Selection(hs[i]).get_selected_particle_indexes():
    #         IMP.atom.Atom(m, pid).set_occupancy(occs[i])

    # IMP.atom.write_multimodel_pdb(hs, str(Path(Path.home(), "xray/tmp/tmp.pdb")))