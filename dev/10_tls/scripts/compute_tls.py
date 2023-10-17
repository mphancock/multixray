from pathlib import Path
import sys

from mmtbx.tls import tools
import mmtbx.model
import iotbx.pdb
from cctbx.array_family import flex

sys.path.append(str(Path(Path.home(), "xray/src")))
import miller_ops

if __name__ == "__main__":
    ###> Get started from PDB
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhf.pdb")
    f_obs_file = Path(Path.home(), "xray/data/reflections/7mhf/7mhf.cif")

    model = mmtbx.model.manager(model_input=iotbx.pdb.input(str(pdb_file)))

    f_obs = miller_ops.get_miller_array(
            f_obs_file=f_obs_file,
            label="_refln.F_meas_au"
        )
    f_obs = miller_ops.clean_miller_array(f_obs)
    status_array = miller_ops.get_miller_array(
        f_obs_file=f_obs_file,
        label="_refln.status"
    )
    flags = status_array.customized_copy(data=status_array.data()=="f")
    f_obs, flags = f_obs.common_sets(other=flags)

    crystal_symmetry = f_obs.crystal_symmetry()
    # xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
    #     crystal_symmetry=crystal_symmetry
    # )

    model.setup_scattering_dictionaries(scattering_table="wk1995")
    model.get_xray_structure().convert_to_isotropic()
    u_iso_start = model.get_xray_structure().extract_u_iso_or_u_equiv()
    # model.get_xray_structure().convert_to_anisotropic()

    # model.process(make_restraints=True)
    # model.setup_scattering_dictionaries(scattering_table="wk1995")
    # model.get_xray_structure().convert_to_isotropic()
    # u_iso_start = model.get_xray_structure().extract_u_iso_or_u_equiv()
    # model.get_xray_structure().convert_to_anisotropic()

    # bounds = [(1,43),(44,100),(101,155),(156,214),(215,292),(293,306)]
    bounds = [(1,100)]
    selections = list()
    for lower_bound, upper_bound in bounds:
        selections.append(model.get_xray_structure().heavy_selection() & model.get_xray_structure().by_index_selection(list(range(lower_bound, upper_bound+1))))

    # selections = [model.get_xray_structure().heavy_selection() & model.get_xray_structure().by_index_selection(list(range(4950)))]
    # selection_strings = ["chain A"]
    # for string in selection_strings:
    #     selections.append(model.selection(string = string))

    # for i in range(selections[0].size()):
    #     print(selections[0][i])

    ################
    selection = flex.bool(model.get_number_of_atoms(), True)

    # selections = []
    # selection_strings = ["chain A", "chain B", "chain C"]
    # for string in selection_strings:
    #     selections.append(model.selection(string = string))
    # ################
    # selection = flex.bool(model.get_number_of_atoms(), True)

    class refinement_flags: pass
    refinement_flags.adp_tls = selections
    model.set_refinement_flags(refinement_flags)
    model.determine_tls_groups(selection_strings=selections, generate_tlsos=selections)
    model.set_refinement_flags(refinement_flags)
    xray_structure = model.get_xray_structure()

    # class refinement_flags: pass
    # refinement_flags.adp_tls = selections
    # model.set_refinement_flags(refinement_flags)
    # model.determine_tls_groups(selection_strings=selections, generate_tlsos=selections)
    # model.set_refinement_flags(refinement_flags)
    # xray_structure = model.get_xray_structure()
    # ################

    # ###> Get TLS <-> Ucart

    T_initial = []
    L_initial = []
    S_initial = []
    T_initial.append([0.11,0.22,0.33,0.12,0.13,0.23])
    L_initial.append([1.11,1.22,1.33,1.12,1.13,1.23])
    S_initial.append([0.11,0.12,0.13,0.21,0.22,0.23,0.31,0.32,-0.33])

    tlsosA = tools.generate_tlsos(selections     = selections,
                                    xray_structure = xray_structure,
                                    T              = T_initial,
                                    L              = L_initial,
                                    S              = S_initial)

    tlsos = tools.generate_tlsos(selections     = selections,
                                xray_structure = xray_structure,
                                T              = T_initial,
                                L              = L_initial,
                                S              = S_initial)

    tlsos = tools.make_tlso_compatible_with_u_positive_definite(
                    tlsos                                       = tlsos,
                    xray_structure                              = xray_structure.deep_copy_scatterers(),
                    selections                                  = selections,
                    max_iterations                              = 50,
                    number_of_u_nonpositive_definite            = 0,
                    eps                                         = 1e-6,
                    number_of_macro_cycles_for_tls_from_uanisos = 30)

    # u_cart_answer = tools.u_cart_from_tls(sites_cart = xray_structure.sites_cart(),
    #                                         selections = selections,
    #                                         tlsos      = tlsos)
    # xray_structure.scatterers().set_u_cart(xray_structure.unit_cell(),
    #                                                                 u_cart_answer)
    tools.show_tls(tlsos = tlsos, text = "ANSWER")

    # T_initial = []
    # L_initial = []
    # S_initial = []
    # T_initial.append([0.11,0.22,0.33,0.12,0.13,0.23])
    # L_initial.append([1.11,1.22,1.33,1.12,1.13,1.23])
    # S_initial.append([0.11,0.12,0.13,0.21,0.22,0.23,0.31,0.32,-0.33])

    # T_initial.append([0.22,0.44,0.66,0.24,0.26,0.46])
    # L_initial.append([2.22,2.44,2.66,2.24,2.26,2.46])
    # S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

    # T_initial.append([0.33,0.66,0.99,0.36,0.39,0.69])
    # L_initial.append([2.33,2.66,2.99,2.36,2.39,2.69])
    # S_initial.append([0.22,0.24,0.26,0.42,0.44,0.46,0.62,0.64,-0.66])

    # tlsosA = tools.generate_tlsos(selections     = selections,
    #                                 xray_structure = xray_structure,
    #                                 T              = T_initial,
    #                                 L              = L_initial,
    #                                 S              = S_initial)

    # tlsos = tools.generate_tlsos(selections     = selections,
    #                             xray_structure = xray_structure,
    #                             T              = T_initial,
    #                             L              = L_initial,
    #                             S              = S_initial)
    # tlsos = tools.make_tlso_compatible_with_u_positive_definite(
    #                 tlsos                                       = tlsos,
    #                 xray_structure                              = xray_structure.deep_copy_scatterers(),
    #                 selections                                  = selections,
    #                 max_iterations                              = 50,
    #                 number_of_u_nonpositive_definite            = 0,
    #                 eps                                         = eps,
    #                 number_of_macro_cycles_for_tls_from_uanisos = 30)

    # u_cart_answer = tools.u_cart_from_tls(sites_cart = xray_structure.sites_cart(),
    #                                         selections = selections,
    #                                         tlsos      = tlsos)
    # xray_structure.scatterers().set_u_cart(xray_structure.unit_cell(),
    #                                                                 u_cart_answer)

    # assert approx_equal(u_cart_answer,
    #         xray_structure.scatterers().extract_u_cart(xray_structure.unit_cell()))


    # tools.show_tls(tlsos = tlsos, text = "ANSWER")

    ###> Set up fmodel
    sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    sfg_params.algorithm = "direct"
    sfg_params.cos_sin_table = False

    fmodel = mmtbx.f_model.manager(xray_structure    = xray_structure,
                                    f_obs             = f_obs,
                                    r_free_flags      = flags,
                                    target_name       = "ls_wunit_k1",
                                    sf_and_grads_accuracy_params = sfg_params)
    fmodel.info(free_reflections_per_bin=250, max_number_of_bins=30).show_all()
    # xray_structure.convert_to_isotropic()
    # xray_structure.set_b_iso(value = 25.0)
    fmodel.update_xray_structure(xray_structure = xray_structure,
                                update_f_calc  = True)
    fmodel.info(free_reflections_per_bin=250, max_number_of_bins=30).show_all()
    print("*"*80)
    ###> TLS refinement against xray data

    for start_tls_value in [None]:#[0.0, tlsosA, None]:
    #for start_tls_value in [None]:
        print(" \n "+str(start_tls_value) + " \n ")
        fmodel_cp = fmodel.deep_copy()
        #for sc in fmodel_cp.xray_structure.scatterers():
        #  sc.flags.set_use_u_aniso(True)
        fmodel_cp.xray_structure.convert_to_anisotropic()

        model.set_xray_structure(fmodel_cp.xray_structure)
        tls_refinement_manager = tools.tls_refinement(
                        fmodel                      = fmodel_cp,
                        model                       = model,
                        selections                  = selections,
                        selections_1d               = None,
                        refine_T                    = 1,
                        refine_L                    = 1,
                        refine_S                    = 1,
                        number_of_macro_cycles      = 1,
                        max_number_of_iterations    = 3,
                        start_tls_value             = start_tls_value,
                        run_finite_differences_test = True,
                        eps                         = 1e-6)
        u_cart = tls_refinement_manager.fmodel.xray_structure.scatterers().extract_u_cart(
                                                            xray_structure.unit_cell())

        for i in range(len(u_cart)):
             print(u_cart[i])

        # if("--comprehensive" in sys.argv[1:]):
        #     format   = "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f"
        #     counter = 0
        #     if(start_tls_value == tlsosA): tolerance = 1.e-6
        #     else: tolerance = 0.02
        #     for m1,m2 in zip(u_cart_answer, u_cart):
        #         counter += 1
        #         if(counter < 10):
        #             print("1=" + format % (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5]))
        #             print("2=" + format % (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]))
        #         assert approx_equal(m1,m2, tolerance)

    # ###> Set up fmodel
    # sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    # sfg_params.algorithm = "direct"
    # sfg_params.cos_sin_table = False
    # dummy = xray_structure.structure_factors(algorithm = sfg_params.algorithm,
    #                                         d_min     = 2.0).f_calc()
    # f_obs = abs(dummy.structure_factors_from_scatterers(
    #                         xray_structure = xray_structure,
    #                         algorithm      = sfg_params.algorithm,
    #                         cos_sin_table  = sfg_params.cos_sin_table).f_calc())
    # flags = f_obs.generate_r_free_flags(fraction=0.01, max_free=2000)

    # fmodel = mmtbx.f_model.manager(xray_structure    = xray_structure,
    #                                 f_obs             = f_obs,
    #                                 r_free_flags      = flags,
    #                                 target_name       = "ls_wunit_k1",
    #                                 sf_and_grads_accuracy_params = sfg_params)
    # fmodel.info(free_reflections_per_bin=250, max_number_of_bins=30).show_all()
    # xray_structure.convert_to_isotropic()
    # xray_structure.set_b_iso(value = 25.0)
    # fmodel.update_xray_structure(xray_structure = xray_structure,
    #                             update_f_calc  = True)
    # fmodel.info(free_reflections_per_bin=250, max_number_of_bins=30).show_all()
    # print("*"*80)
    # ###> TLS refinement against xray data
    # if (not "--comprehensive" in sys.argv[1:]):
    #         number_of_macro_cycles   = 1
    #         max_number_of_iterations = 3
    # else:
    #         number_of_macro_cycles   = 100
    #         max_number_of_iterations = 50

    # for start_tls_value in [None]:#[0.0, tlsosA, None]:
    # #for start_tls_value in [None]:
    #     print(" \n "+str(start_tls_value) + " \n ")
    #     fmodel_cp = fmodel.deep_copy()
    #     #for sc in fmodel_cp.xray_structure.scatterers():
    #     #  sc.flags.set_use_u_aniso(True)
    #     fmodel_cp.xray_structure.convert_to_anisotropic()

    #     if(start_tls_value is None):
    #         run_finite_differences_test = True
    #     else: run_finite_differences_test = False
    #     model.set_xray_structure(fmodel_cp.xray_structure)
    #     tls_refinement_manager = tools.tls_refinement(
    #                     fmodel                      = fmodel_cp,
    #                     model                       = model,
    #                     selections                  = selections,
    #                     selections_1d               = None,
    #                     refine_T                    = 1,
    #                     refine_L                    = 1,
    #                     refine_S                    = 1,
    #                     number_of_macro_cycles      = number_of_macro_cycles,
    #                     max_number_of_iterations    = max_number_of_iterations,
    #                     start_tls_value             = start_tls_value,
    #                     run_finite_differences_test = run_finite_differences_test,
    #                     eps                         = eps)
    #     u_cart = tls_refinement_manager.fmodel.xray_structure.scatterers().extract_u_cart(
    #                                                         xray_structure.unit_cell())
    #     if("--comprehensive" in sys.argv[1:]):
    #         format   = "%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f"
    #         counter = 0
    #         if(start_tls_value == tlsosA): tolerance = 1.e-6
    #         else: tolerance = 0.02
    #         for m1,m2 in zip(u_cart_answer, u_cart):
    #             counter += 1
    #             if(counter < 10):
    #                 print("1=" + format % (m1[0],m1[1],m1[2],m1[3],m1[4],m1[5]))
    #                 print("2=" + format % (m2[0],m2[1],m2[2],m2[3],m2[4],m2[5]))
    #             assert approx_equal(m1,m2, tolerance)