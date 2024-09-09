from pathlib import Path
import sys
import argparse
import numpy as np
import pandas as pd
import random
import math
import shutil

import IMP
import IMP.atom
import IMP.core
import IMP.isd
import IMP.algebra

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import xray_restraint
import trackers
import com_restraint
import log_statistics
import align_imp
import pdb_writer
import com_optimizer_state
import update_weights_optimizer_state
import reset
import molecular_dynamics
from params import write_params_txt, write_params_csv, read_job_csv
import weight_restraint
import miller_ops
import weights
import utility
import multi_state_multi_condition_model
from derivatives import evaluate_df_dict


if __name__ == "__main__":
    m = IMP.Model()
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhf_refine.pdb")
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()

    ### SCORING
    rs = list()
    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmm.charmm_restraints(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)

    # Is turning this off going to be really bad?
    # rset_charmm.set_weight(1/N)
    rs.append(rset_charmm)

    test_pid = pids[0]
    test_xyz = IMP.core.XYZ(m, test_pid)

    print("dxyz: ", test_xyz.get_derivatives())

    print("EVALUATING CHARMM RESTRAINT SET")
    rset_charmm.evaluate(calc_derivs=True)

    print("dxyz: ", test_xyz.get_derivatives())
    # print("xyz: ", test_xyz.get_coordinates())



    md = IMP.atom.MolecularDynamics(m)
    ps_0 = [m.get_particle(pid) for pid in pids]
    # s_v = IMP.atom.VelocityScalingOptimizerState(m, ps_0, T_0)
    # md.add_optimizer_state(s_v)
    md.set_particles(ps_0)
    sf = IMP.core.RestraintsScoringFunction([rset_charmm])
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    md.set_temperature(300)
    md.setup(ps_0)
    for pid in pids:
        IMP.atom.LinearVelocity.setup_particle(m, pid)

    md.assign_velocities(300)
    md.set_maximum_time_step(1)

    print("UPDATING")
    md.simulate(1)

    print("dxyz: ", test_xyz.get_derivatives())
    # print("xyz: ", test_xyz.get_coordinates())

    print("CHECKING DERIVATIVES")
    # rset_charmm.evaluate(True)
    # rset_charmm.evaluate(True)

    # print("dxyz: ", test_xyz.get_derivatives())
    # print("xyz: ", test_xyz.get_coordinates())


    df_dict = evaluate_df_dict(
        m=m,
        pids=pids,
        r=rset_charmm
    )
    print("dxyz: ", df_dict[test_pid])
    print("dxyz: ", test_xyz.get_derivatives())
