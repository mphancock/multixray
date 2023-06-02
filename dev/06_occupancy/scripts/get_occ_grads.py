from pathlib import Path
import sys
import numpy as np

import IMP
import IMP.atom

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import iotbx

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
import miller_ops
import xray_struct


def get_r_factor(
        xray_structure
):
    cif_file = Path(Path.home(), "xray/data/reflections/7mhk/7mhk.cif")
    f_obs_array = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.F_meas_au"
    )
    f_obs_array = miller_ops.clean_miller_array(f_obs_array)

    # Set flags.
    status_array = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.status"
    )
    flags_array = status_array.customized_copy(data=status_array.data()=="f")
    f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

    # Compute r_factor
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

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs_array,
        r_free_flags=flags_array,
        target_name="ml"
    )
    f_model_manager.update_all_scales()

    fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)
    fmodels.update_xray_structure(
        xray_structure=xray_structure,
        update_f_calc=True)

    r_work = f_model_manager.r_work()
    r_free = f_model_manager.r_free()
    r_all = f_model_manager.r_all()

    fmodels_target_and_gradients = fmodels.target_and_gradients(
        compute_gradients=True)
    score = fmodels_target_and_gradients.target()
    grads = fmodels_target_and_gradients.gradients()

    occ_grads = list()
    for i in range(len(xray_structure.scatterers())):
        occ_id = 3+i*4
        occ_grads.append(grads[occ_id])

    print(len(occ_grads))
    print(occ_grads[0])
    print(occ_grads[1])

    print(len(grads))
    print(np.mean(grads))

    print(r_work, r_free)

if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhk/7mhk_clean_h20.pdb")
    uc = "114.300, 54.290, 44.970, 90.00, 102.12, 90.00"
    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=uc,
        space_group_symbol="C 1 2 1"
    )

    xray_structure_cctbx = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    n_state = 1
    m = IMP.Model()
    hs = list()
    for i in range(n_state):
        h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        hs.append(h)

    xray_structure_imp = xray_struct.get_xray_structure(
        m,
        uc_dim=uc,
        sg_symbol="C 1 2 1"
    )

    get_r_factor(xray_structure_imp)
