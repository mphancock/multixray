from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import argparse
import time
import shutil

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import params
import trackers
import log_statistics
import pdb_writer
from utility import pool_read_pdb
from derivatives import DerivScoreState
sys.path.append(str(Path(Path.home(), "xray/src")))
from params import write_params_txt, write_params_csv, read_job_csv
from miller_ops import get_crystal_symmetry, get_f_obs, get_flags, clean_miller_array
import multi_state_multi_condition_model
from charmm import get_charmm_restraint_set, CHARMMDerivHolder
import xray_restraint


if __name__ == "__main__":
    pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/271_native_2_wxray/1/output_5/pdbs/291.pdb")
    ref_pdb_file = Path(Path.home(), "xray/dev/45_synthetic_native_4/data/pdbs/native.pdb")
    w_mat = np.array([[0.75, 0.25], [0.25, 0.75]])
    ref_w_mat = np.array([[0.75, 0.25], [0.25, 0.75]])
    cif_files = [Path(Path.home(), "xray/dev/45_synthetic_native_4/data/cifs/native_2_0.cif"), Path(Path.home(), "xray/dev/45_synthetic_native_4/data/cifs/native_2_1.cif")]

    N, J = w_mat.shape

    crystal_symmetries = []
    for cond in range(J):
        crystal_symmetries.append(get_crystal_symmetry(
            cif_file=cif_files[cond]
        ))

    msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
        pdb_files=[pdb_file],
        w_mat=w_mat,
        crystal_symmetries=crystal_symmetries
    )

    ### SCORING
    m, hs = msmc_m.get_m(), msmc_m.get_hs()
    r_charmm = get_charmm_restraint_set(m, hs)

    # cif file here is a string.
    rset_xray = IMP.RestraintSet(m, 1.0)
    charmm_deriv_holder = CHARMMDerivHolder()
    for i in range(len(cif_files)):
        cif_file = cif_files[i]
        if not cif_file:
            continue

        f_obs = get_f_obs(cif_file=cif_file)
        f_obs = clean_miller_array(f_obs)

        flags = get_flags(cif_file=cif_file)
        f_obs, flags = f_obs.common_sets(other=flags)

        r_xray = xray_restraint.XtalRestraint(
            msmc_m=msmc_m,
            cond=i,
            f_obs=f_obs,
            free_flags=flags,
            w_xray=1,
            update_scale=True,
            update_k1=True,
            update_freq=1,
            charmm_holder=charmm_deriv_holder,
            ref_com=None
        )

        rset_xray.add_restraint(r_xray)

    r_xrays = [rset_xray.get_restraint(i) for i in range(rset_xray.get_number_of_restraints())]

    deriv_score_state = DerivScoreState(
        m=m,
        pids=msmc_m.get_pids(),
        charmm_holder=charmm_deriv_holder,
        xray_rs=r_xrays,
        w_xray=1
    )
    m.add_score_state(deriv_score_state)

    for r_xray in r_xrays:
        r_xray.evaluate(True)

    sf_refine = IMP.core.RestraintsScoringFunction([r_charmm, rset_xray])
    cg = IMP.core.ConjugateGradients(msmc_m.get_m())
    cg.set_scoring_function(sf_refine)
    cg.optimize(50)

    print("r_free", r_xrays[0].get_r_free())
    print("r_free", r_xrays[1].get_r_free())
    ff = r_charmm.evaluate(True)
    print("ff", ff)