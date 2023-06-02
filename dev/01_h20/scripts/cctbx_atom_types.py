from pathlib import Path
import sys

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
import cctbx_score


def get_r_factor(
        xray_structure
):
    cif_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")
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

    print(r_work, r_free)

if __name__ == "__main__":
    # pdb_file = Path(Path.home(), "xray/dev/01_h20/data/3ca7_only_A.pdb")
    # pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7.pdb")
    # pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean_h20.pdb")
    # pdb_file = Path(Path.home(), "xray/dev/01_h20/data/3ca7_test.pdb")

    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell="58.3050, 36.1540, 25.3620, 90.0000, 103.0900, 90.0000",
        space_group_symbol="C 1 2 1"
    )

    xray_structure_cctbx = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    # xray_structure_imp = xray_struct.get_cctbx_multi_structure_from_pdbs(
    #     pdb_files=[pdb_file],
    #     weights=[1],
    #     crystal_symmetry=crystal_symmetry
    # )

    n_state = 1
    m = IMP.Model()
    hs = list()
    for i in range(n_state):
        h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        hs.append(h)

    xray_structure_imp = xray_struct.get_xray_structure(
        m,
        uc_dim="58.3050, 36.1540, 25.3620, 90.0000, 103.0900, 90.0000",
        sg_symbol="C 1 2 1"
    )

    # xray_structure_cctbx.show_scatterer_flags_summary()
    # xray_structure_imp.show_scatterer_flags_summary()

    for i in range(xray_structure_cctbx.scatterers().size()):
        scatterer_cctbx = xray_structure_cctbx.scatterers()[i]
        scatterer_imp = xray_structure_imp.scatterers()[i]

        if scatterer_cctbx.label != scatterer_imp.label:
            print("ERROR: labels do not match.")
            print(scatterer_cctbx.label, scatterer_imp.label)

        if scatterer_cctbx.site != scatterer_imp.site:
            print("ERROR: sites do not match.")
            print(scatterer_cctbx.site, scatterer_imp.site)

        if scatterer_cctbx.occupancy != scatterer_imp.occupancy:
            print(scatterer_cctbx.label, scatterer_imp.label)
            print("ERROR: occupancies do not match.")
            print(scatterer_cctbx.occupancy, scatterer_imp.occupancy)

        if scatterer_cctbx.scattering_type != scatterer_imp.scattering_type:
            print("ERROR: scattering types do not match.")
            print(scatterer_cctbx.scattering_type, scatterer_imp.scattering_type)

        if scatterer_cctbx.fp != scatterer_imp.fp:
            print("ERROR: f' do not match.")
            print(scatterer_cctbx.fp, scatterer_imp.fp)

        if scatterer_cctbx.fdp != scatterer_imp.fdp:
            print("ERROR: f'' do not match.")
            print(scatterer_cctbx.fdp, scatterer_imp.fdp)

        if scatterer_cctbx.u_iso != scatterer_imp.u_iso:
            print("ERROR: u_iso do not match.")
            print(scatterer_cctbx.u_iso, scatterer_imp.u_iso)

        if scatterer_cctbx.u_star != scatterer_imp.u_star:
            print("ERROR: u_star do not match.")
            print(scatterer_cctbx.u_star, scatterer_imp.u_star)

        # if scatterer_cctbx.flags != scatterer_imp.flags:
        #     print("ERROR: flags do not match.")
        #     print(scatterer_cctbx.flags, scatterer_imp.flags)

        if scatterer_cctbx.multiplicity() != scatterer_imp.multiplicity():
            print("ERROR: multiplicity do not match.")
            print(scatterer_cctbx.multiplicity(), scatterer_imp.multiplicity())

        if scatterer_cctbx.weight_without_occupancy() != scatterer_imp.weight_without_occupancy():
            print("ERROR: weight_without_occupancy do not match.")
            print(scatterer_cctbx.weight_without_occupancy(), scatterer_imp.weight_without_occupancy())

    get_r_factor(xray_structure_cctbx)
    get_r_factor(xray_structure_imp)
