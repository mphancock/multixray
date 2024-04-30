from pathlib import Path
import pickle

import mmtbx.f_model
import mmtbx.model
import cctbx.crystal
import cctbx.xray
import iotbx

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
sys.path.append(str(Path(Path.home(), "xray/data/cifs/scripts")))
# sys.path.append("../src")
import miller_ops
import generate_fmodel


if __name__ == "__main__":
    score_dir = Path(Path.home(), "xray/dev/32_score")

    pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/183_test/788442/output_17/pdbs/123.pdb")
    occs = [0.5429301174113870, 0.4570698825886140]
    # pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/54_7mhf_decoys_100/4545474/output_3701/pdbs/49.pdb")
    # occs = [0.5099659036431840, 0.4900340963568160]
    # pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/54_7mhf_decoys_100/4545474/output_4676/pdbs/1.pdb")
    # occs = [0.4885097353792970, 0.511490264620703]
    # pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/54_7mhf_decoys_100/4545474/output_4564/pdbs/23.pdb")
    # occs = [0.7395988883650760, 0.260401111634924]

    f_obs_file = Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/0/0.cif")

    # pdb_file = Path(score_dir, "data/decoy.pdb")
    # f_obs_file = Path(score_dir, "data/0.cif")

    # pdb_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs/2_state_0/0.pdb")
    # f_obs_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/2_state_0_noise/0.cif")

    f_obs = miller_ops.get_miller_array(
            f_obs_file=f_obs_file,
            label="_refln.F_meas_au"
        )

    print(f_obs.indices()[0], f_obs.data()[0])

    f_obs = miller_ops.clean_miller_array(f_obs)
    status_array = miller_ops.get_miller_array(
        f_obs_file=f_obs_file,
        label="_refln.status"
    )
    flags = status_array.customized_copy(data=status_array.data()=="f")
    f_obs, flags = f_obs.common_sets(other=flags)

    d_min = 2.0
    f_obs = miller_ops.filter_f_obs_resolution(
            f_obs=f_obs,
            d_max=None,
            d_min=d_min
    )
    flags = miller_ops.filter_f_obs_resolution(
            f_obs=flags,
            d_max=None,
            d_min=d_min
    )

    print(f_obs.indices()[0], f_obs.data()[0])


    crystal_symmetry = f_obs.crystal_symmetry()
    model = mmtbx.model.manager(
        model_input=iotbx.pdb.input(str(pdb_file)),
        crystal_symmetry=crystal_symmetry
    )

    xray_structure = model.get_xray_structure()


    n_scatt = int(xray_structure.scatterers().size())
    print(type(n_scatt))
    for i in range(n_scatt//2):
        xray_structure.scatterers()[i].occupancy = occs[0]
        xray_structure.scatterers()[i+n_scatt//2].occupancy = occs[1]

    xray_structure.scatterers().flags_set_grads(
        state=False
    )
    xray_structure.scatterers().flags_set_grad_site(
        iselection=xray_structure.all_selection().iselection()
    )
    xray_structure.scatterers().flags_set_grad_occupancy(
        iselection=xray_structure.all_selection().iselection()
    )

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        r_free_flags=flags,
        target_name="ml",
        max_number_of_bins=1
    )
    f_model_manager.update_all_scales(apply_scale_k1_to_f_obs=False,remove_outliers=False)

    f_obs = f_model_manager.f_obs()
    # print(len(f_obs.indices()), len(f_obs.data()))
    # print(f_obs.indices()[0], f_obs.data()[0])

    for i in range(len(f_obs.data())):
        print(f_model_manager.f_obs().indices()[i],f_model_manager.f_calc().amplitudes().data()[i],f_model_manager.f_model().amplitudes().data()[i], f_model_manager.f_obs().data()[i])

        break
    # print(f_model_manager.f_model_no_scales().data()[-1],f_model_manager.f_model().data()[-1], f_model_manager.f_obs().data()[-1])

    # print(f_model_manager.k_masks()[0][0])
    # print(len(f_model_manager.k_isotropic()), f_model_manager.k_isotropic()[0],f_model_manager.k_isotropic()[-1])
    # print(f_model_manager.k_anisotropic()[0],f_model_manager.k_anisotropic()[-1])
    # print(f_model_manager.k_sol)
    # print(f_model_manager.b_sol)
    # print(f_model_manager.b_cart)

    # f_model_manager.show_short()
    f_model_manager.show()


    # f = open('alpha_beta', 'rb')
    # alpha_beta = pickle.load(f)

    # beta = alpha_beta[1]
    # for i in range(len(beta.data())):
    #     beta.data()[i] = beta.data()[i] * 1000

    # f_model_manager.alpha_beta_cache = alpha_beta


    alpha_beta = f_model_manager.alpha_beta()
    print(f_model_manager.alpha_beta()[0].size(), f_model_manager.f_obs().size())
    print(alpha_beta[0].data()[0], alpha_beta[1].data()[0])
    print(alpha_beta[0].data()[-1], alpha_beta[1].data()[-1])

    f = open('/wynton/home/sali/mhancock/xray/dev/32_score/alpha_beta', 'ab')
    b = pickle.dump(alpha_beta, f)
    f.close()

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

    # print(xray_structure.scatterers()[0].occupancy, xray_structure.scatterers()[-1].occupancy)
    print(r_free, r_work)
    print("likelihood:", score)