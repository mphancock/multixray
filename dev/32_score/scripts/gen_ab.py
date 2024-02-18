from pathlib import Path
import pickle
import pandas as pd

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
    cif_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/7mhf/0")
    native_df = pd.read_csv(Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/7mhf.csv"), index_col=0)

    for i in range(len(native_df)):
        pdb_file = native_df.loc[i, "pdb"]
        ab_file = Path(score_dir, "data/ab/ab_{}".format(i))
        ab_file.unlink(missing_ok=True)

        f_obs_file = Path(cif_dir, "{}.cif".format(i))

        occ_0 = native_df.loc[i, "weight_0_0"]
        occ_1 = native_df.loc[i, "weight_0_1"]

        f_obs = miller_ops.get_miller_array(
                f_obs_file=f_obs_file,
                label="_refln.F_meas_au"
            )

        n_obs = len(f_obs.data())
        print(n_obs)

        f_obs = miller_ops.clean_miller_array(f_obs)
        status_array = miller_ops.get_miller_array(
            f_obs_file=f_obs_file,
            label="_refln.status"
        )

        # for i in range(10):
        #     print(f_obs.indices()[i], f_obs.data()[i], status_array.data()[i])

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

        # print(xray_structure.show_scatterers())
        n_scatt = int(xray_structure.scatterers().size())
        print(type(n_scatt))
        for i in range(n_scatt//2):
            xray_structure.scatterers()[i].occupancy = occ_0
            xray_structure.scatterers()[i+n_scatt//2].occupancy = occ_1

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
        # for i in range(len(f_obs.data())):
        #     print(f_model_manager.f_obs().indices()[i],f_model_manager.f_calc().amplitudes().data()[i],f_model_manager.f_model().amplitudes().data()[i], f_model_manager.f_obs().data()[i], status_array.data()[i], flags.data()[i])

        f_model_manager.show()

        alpha_beta = f_model_manager.alpha_beta()
        print(f_model_manager.alpha_beta()[0].size(), f_model_manager.f_obs().size())
        print(alpha_beta[0].data()[0], alpha_beta[1].data()[0])
        print(alpha_beta[0].data()[-1], alpha_beta[1].data()[-1])

        f = open(ab_file, 'ab')
        b = pickle.dump(alpha_beta, f)
        f.close()

        r_work = f_model_manager.r_work()
        r_free = f_model_manager.r_free()
        r_all = f_model_manager.r_all()
        print(r_free, r_work)

        # fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)

        # fmodels_target_and_gradients = fmodels.target_and_gradients(compute_gradients=True)
        # score = fmodels_target_and_gradients.target()
        # grads = fmodels_target_and_gradients.gradients()
        # print("likelihood:", score)

        break