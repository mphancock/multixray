import sys
from pathlib import Path
import shutil
import numpy as np
import math

import IMP
import IMP.core
import IMP.atom

sys.path.append("../../../src")
import charmm
import params
import trackers
import log_statistics
import pdb_writer
from utility import pool_read_pdb
from multi_state_multi_condition_model import MultiStateMultiConditionModel


if __name__ == "__main__":
    msmc_model = MultiStateMultiConditionModel(
        pdb_file=Path(Path.home(), "Documents/xray/data/pdbs/3k0m/3k0n.pdb"),
        w_mat=np.array([[1]])
    )

    hs = msmc_model.get_hs()
    m = msmc_model.get_m()

    pid = IMP.atom.Selection(hs[0]).get_selected_particles()[0]
    d = IMP.core.XYZR(m, pid)
    print(d.get_coordinates())

    for h in hs:
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

        # If the initial structure is non-physiological, don't refine
        ff_cur = rset_charmm.evaluate(False)

        o_states = list()

        # if log_file:
        #     log_ostate = log_statistics.LogStatistics(
        #         m=m,
        #         all_trackers=all_trackers,
        #         log_file=log_file,
        #         log_freq=1000
        #     )
        #     o_states.append(log_ostate)

        sf = IMP.core.RestraintsScoringFunction(rs)
        cg = IMP.core.ConjugateGradients(m)
        cg.set_scoring_function(sf)

        for o_state in o_states:
            cg.add_optimizer_state(o_state)

        cg.optimize(50)

        ff_new = rset_charmm.evaluate(False)
        print(ff_cur, ff_new)

    print(d.get_coordinates())

    msmc_model.write_pdb_file(Path(Path.home(), "Documents/xray/data/pdbs/3k0m/3k0n_refine.pdb"))
