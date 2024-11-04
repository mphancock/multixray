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

from cctbx.crystal import symmetry

sys.path.append(str(Path(Path.home(), "xray/src")))
from charmm import get_charmm_restraint_set, CHARMMDerivHolder
import xray_restraint
import trackers
import com_restraint
import log_statistics
import align_imp
import pdb_writer
import com_optimizer_state
import update_weights_optimizer_state
import reset
from simulated_annealing import SimulatedAnnealing, SimulatedAnnealingSchedule
from params import write_params_txt, write_params_csv, read_job_csv
from miller_ops import get_crystal_symmetry, get_f_obs, get_flags, clean_miller_array
import weights
import utility
import multi_state_multi_condition_model
from derivatives import evaluate_df_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_file")
    parser.add_argument("--n_pdb", type=int)
    args = parser.parse_args()


    # pdb_files = list()
    # for job_id in range(36):
    #     job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/263_sb_sa/{}".format(job_id))
    #     for out_dir in job_dir.glob("output*"):
    #         pdb_file = Path(out_dir, "pdbs/500.pdb")

    #         if pdb_file.exists():
    #             pdb_files.append(pdb_file)

    # import multiprocessing

    # pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    # pool_results = pool_obj.imap(refine, pdb_files)

    # r_free_df = pd.DataFrame()
    # for pool_result in pool_results:
    #     score_dict = pool_result

    #     pdb_file = score_dict["pdb_file"]

    #     job_id = int(pdb_file.parts[-4])
    #     out_id = int(pdb_file.parts[-3].split("_")[-1])

    #     r_free_df.loc[job_id, out_id] = score_dict["r_free_0"]

    # r_free_df.to_csv("../data/r_free_100.csv")


        # # pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/263_sb_sa/2/output_0/pdbs/100.pdb")
        #     charmm, r_work_0, refine(pdb_file)