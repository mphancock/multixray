import sys
from pathlib import Path
import shutil

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))

import charmm
import params
import trackers
import log_statistics
import pdb_writer


def refine(
        hs,
        n_step
):
    m = hs[0].get_model()
    for i in range(len(hs)):
        h = hs[i]
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        print(len(pids))
        ps = [m.get_particle(pid) for pid in pids]

        rs = list()
        rset_charmm = IMP.RestraintSet(m, 1.0)
        charmm_rs = charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
        rset_charmm.add_restraints(charmm_rs)
        rset_charmm.set_weight(1)
        rs.append(rset_charmm)

        all_trackers = list()
        step_tracker = trackers.StepTracker(
            name="step",
            m=m
        )
        all_trackers.append(step_tracker)

        time_tracker = trackers.TimeTracker(
            name="time",
            m=m
        )
        all_trackers.append(time_tracker)

        ff_tracker = trackers.fTracker(
            name="ff",
            r=rset_charmm
        )
        all_trackers.append(ff_tracker)

        o_states = list()

        sf = IMP.core.RestraintsScoringFunction(rs)
        cg = IMP.core.ConjugateGradients(m)
        cg.set_scoring_function(sf)

        for o_state in o_states:
            cg.add_optimizer_state(o_state)

        cg.optimize(n_step)

    return 0

if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/tmp/9.pdb")
    print(pdb_file)

    out_pdb_file = Path(Path.home(), "xray/tmp/9_ref.pdb")
    n_steps = 100
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    refine(
        hs=hs,
        n_step=n_steps
    )

    IMP.atom.write_multimodel_pdb(hs, str(out_pdb_file))




