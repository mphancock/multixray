from pathlib import Path
import sys
import numpy as np

import IMP
import IMP.atom
import IMP.isd

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import iotbx
import mmtbx.refinement.occupancies

from libtbx.utils import null_out
from mmtbx import monomer_library
import mmtbx.model
import iotbx.pdb

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
import miller_ops
import xray_struct
import weight_restraint

from libtbx.test_utils import show_diff, Exception_expected
from libtbx.phil import interface
import libtbx.load_env
import libtbx.phil
from six.moves import cStringIO as StringIO
import sys

from cctbx.array_family import flex


def set_occs(
    hs,
    occs
):
    for i in range(len(hs)):
        h = hs[i]
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        for pid in pids:
            IMP.atom.Atom(m, pid).set_occupancy(occs[i])


if __name__ == "__main__":
    # print("HELLO")

    # if (not libtbx.env.has_module(name="phenix")):
    #     print("phenix module not available: skipping advanced tests")

    fmodels = None
    selections = None
    constrained_groups_selections = None
    par_initial = None
    max_number_of_iterations = 1


    f_obs_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7_refine_2.cif")
    # Set f_obs.
    f_obs_array = miller_ops.get_miller_array(
        f_obs_file=f_obs_file,
        label="_refln.F_meas_au"
    )
    f_obs_array = miller_ops.clean_miller_array(f_obs_array)

    # Set flags.
    status_array = miller_ops.get_miller_array(
        f_obs_file=f_obs_file,
        label="_refln.status"
    )
    flags_array = status_array.customized_copy(data=status_array.data()=="f")
    f_obs, flags_array = f_obs_array.common_sets(other=flags_array)

    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_refine_2.pdb")

    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    occs = [.85,.15]
    set_occs(
        hs=hs,
        occs=occs
    )

    crystal_symmetry = f_obs.crystal_symmetry()

    ws = list()
    pids_0 = IMP.atom.Selection(hs[0]).get_selected_particle_indexes()[0]
    at = IMP.atom.Atom(m, pids_0)
    w = IMP.isd.Weight.setup_particle(m, pids_0)
    # w.add_weight(occs[0])
    # w.add_weight(occs[1])
    w.set_number_of_weights(2)
    w.set_weights([occs[0],occs[1]])
    w.set_weights_are_optimized(True)

    print(w.get_weight(0), w.get_weight(1))

    wr = weight_restraint.WeightRestraint(
        m=m,
        hs=hs,
        w=w,
        f_obs=f_obs,
        flags=flags_array,
        scale=1.0
    )
    sf = IMP.core.RestraintsScoringFunction([wr])
    # sf.evaluate(True)

    cg = IMP.core.ConjugateGradients(m)
    cg.set_scoring_function(sf)

    rg = IMP._get_range(m, w.get_weight_key(0))
    print("range: ", rg)
    # cg.set_max_change(.2)
    # cg.set_gradient_threshold(1.0)
    cg.optimize(25)
    print(wr.get_best_occs())

    # cg.set_max_change(.1)
    # cg.optimize(10)

    # wr.set_scale(.1)
    # cg.optimize(10)



    # for i in range(100):
    #     xray_structure = xray_struct.get_xray_structure(
    #         m=m,
    #         crystal_symmetry=crystal_symmetry
    #     )

    #     # m = mmtbx.model.manager(
    #     #     model_input = iotbx.pdb.input(file_name=str(pdb_file))
    #     # )
    #     # xray_structure = m.get_xray_structure()
    #     xray_structure.scatterers().flags_set_grads(
    #         state=False
    #     )
    #     xray_structure.scatterers().flags_set_grad_site(
    #         iselection=xray_structure.all_selection().iselection()
    #     )
    #     xray_structure.scatterers().flags_set_grad_occupancy(
    #         iselection=xray_structure.all_selection().iselection()
    #     )

    #     # r_factor_target = False

    #     f_model_manager = mmtbx.f_model.manager(
    #         xray_structure=xray_structure,
    #         f_obs=f_obs,
    #         r_free_flags=flags_array,
    #         target_name="ml"
    #     )
    #     f_model_manager.update_all_scales()
    #     fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)

    #     fmodels.create_target_functors()
    #     fmodels.prepare_target_functors_for_minimization()

    #     fmodels.update_xray_structure(update_f_calc = True)
    #     fmodels_target_and_gradients = fmodels.target_and_gradients(
    #     weights=None,
    #     compute_gradients=True,
    #     occupancy=True
    #     )

    #     f = fmodels_target_and_gradients.target()
    #     g =  fmodels_target_and_gradients.gradients()

    #     print(g.size())

    #     occ_gs = list()
    #     for i in range(2):
    #         occ_gs.append(list())
    #         for j in range(402):
    #             # print(g[i*402+j])
    #             occ_gs[i].append(g[i*402+j])

    #     print(occ_gs[0][0])

    #     g0, g1 = sum(occ_gs[0])/402, sum(occ_gs[1])/402
    #     # g0, g1 = sum(occ_gs[0]), sum(occ_gs[1])

    #     g0 = g0 - g1
    #     g1 = 0

    #     if abs(g0) > .1:
    #         g0 = .1

    #     print(g0, g1)

    #     occs = [occs[0]-g0,occs[1]-g1]
    #     occs = [occs[0]/sum(occs),occs[1]/sum(occs)]
    #     print(occs)
    #     print(f)

    #     set_occs(
    #         hs=hs,
    #         occs=occs
    #     )




    #     # fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)

    #     # occ_sel = mmtbx.refinement.occupancies.occupancy_selections(
    #     #     m,
    #     #     other_individual_selection_strings=['resseq 49 or resseq 50']
    #     # )

    #     # for i in range(len(occ_sel)):
    #     #     print("NEW GROUP")
    #     #     print(type(occ_sel[i]))
    #     #     for j in range(len(occ_sel[i])):
    #     #         print("NEW SUBGROUP")
    #     #         print(type(occ_sel[i][j]))
    #     #         for k in range(occ_sel[i][j].size()):
    #     #             print(occ_sel[i][j][k])
    #     #             print(type(occ_sel[i][j][k]))

