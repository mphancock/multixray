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
from utility import pool_read_pdb


"""
I read and refine in a single function to reduce memory usage versus loading all hierarchies/models at once.
"""
def read_pdb_and_refine(
    pool_params
):
    pdb_file = pool_params["pdb_file"]
    out_pdb_file = pool_params["out_pdb_file"]
    n_step = pool_params["n_step"]
    log_file = pool_params["log_file"]

    if pdb_file.exists():
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    else:
        return RuntimeError("{} file does not exist.".format(pdb_file))

    refine_hs(hs=hs, n_step=n_step, log_file=log_file)
    IMP.atom.write_multimodel_pdb(hs, str(out_pdb_file))

    return out_pdb_file


def pool_refine(
    pool_param
):
    hs = pool_param["hs"]
    n_step = pool_param["n_step"]
    log_file = pool_param["log_file"]

    return refine_hs(hs=hs, n_step=n_step, log_file=log_file)


def refine_hs(
        hs,
        n_step,
        log_file
):
    m = hs[0].get_model()
    for i in range(len(hs)):
        h = hs[i]
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        # print(len(pids))
        for pid in pids:
            IMP.core.XYZ(m, pid).set_coordinates_are_optimized(True)

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

        if log_file:
            log_ostate = log_statistics.LogStatistics(
                m=m,
                all_trackers=all_trackers,
                log_file=log_file,
                log_freq=1000
            )
            o_states.append(log_ostate)

        sf = IMP.core.RestraintsScoringFunction(rs)
        cg = IMP.core.ConjugateGradients(m)
        cg.set_scoring_function(sf)

        for o_state in o_states:
            cg.add_optimizer_state(o_state)

        cg.optimize(n_step)

    # print(IMP.core.XYZ(m, IMP.atom.Selection(hs[0]).get_selected_particle_indexes()[0]).get_x())

    return m, hs
