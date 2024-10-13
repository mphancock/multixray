from pathlib import Path

import numpy as np

import IMP
import IMP.atom
import IMP.core

import sys
sys.path.append("..")
from multi_state_multi_condition_model import MultiStateMultiConditionModel
from charmm import get_charmm_restraint_set


if __name__ == "__main__":
    w_mat = np.array([[1.0]])
    # pdb_file = Path("/wynton/home/sali/mhancock/xray/tmp/min.pdb")
    pdb_file = Path("/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_no_dms.pdb")

    msmc_m = MultiStateMultiConditionModel(
        pdb_files=[pdb_file],
        w_mat=w_mat,
        crystal_symmetries=None
    )

    h, m = msmc_m.get_hs()[0], msmc_m.get_m()

    # atoms = IMP.atom.get_by_type(h, IMP.atom.ATOM_TYPE)
    # cont = IMP.container.ListSingletonContainer(m, atoms)
    # nbl = IMP.container.ClosePairContainer(cont, 5, 0)

    # r = CHARMMRestraint(msmc_m=msmc_m)
    # sf = IMP.core.RestraintsScoringFunction([r])

    rset_charmm = get_charmm_restraint_set(m, [h])
    sf = IMP.core.RestraintsScoringFunction([rset_charmm])

    # IMP.set_log_level(IMP.VERBOSE)
    # charmm_rs[-2].evaluate(True)

    # m.update()

    # IMP.set_log_level(IMP.VERBOSE)
    sf.evaluate(True)

    # pid = IMP.atom.Selection(h).get_selected_particle_indexes()[0]

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    ps = [m.get_particle(pid) for pid in pids]

    # charmm_rs[-1].do_add_score_and_derivatives(IMP.ScoreAccumulator())

    cg = IMP.core.ConjugateGradients(m)
    cg.set_scoring_function(sf)
    cg.optimize(1000)

    for r in rset_charmm.get_restraints():
        print(r.get_name(), r.evaluate(True))


    md = IMP.atom.MolecularDynamics(m)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)
    md.setup(ps)
    md.set_temperature(300)
    md.set_maximum_time_step(1.0)
    md.assign_velocities(300)
    vel_therm = IMP.atom.BerendsenThermostatOptimizerState(ps, 300, 10)
    md.add_optimizer_state(vel_therm)

    o_state = IMP.atom.WritePDBOptimizerState(m, ps, "/wynton/home/sali/mhancock/xray/tmp/out/%1%.pdb")
    o_state.set_period(1000)
    md.add_optimizer_state(o_state)
    md.simulate(100000)

    IMP.atom.write_pdb(h, str(Path(Path.home(), "xray/tmp/out.pdb")))

    for r in rset_charmm.get_restraints():
        print(r.get_name(), r.evaluate(True))

    # print(nbl.get_number_of_full_rebuilds())






