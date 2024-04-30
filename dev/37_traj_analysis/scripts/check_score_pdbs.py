from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import os

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
# from refine import refine, pool_refine
from utility import pool_read_pdb


if __name__ == "__main__":
    exp_name = "179_exp"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    analysis_dir = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/7mhf", exp_name)

    score_files = list(analysis_dir.glob("ref_25000*.csv"))

    empty_files = 0
    for score_file in score_files:
        print(score_file)
        score_df = pd.read_csv(score_file, index_col=0)
        contains_nan = score_df['r_free'].isna().any()

        if contains_nan:
            empty_files += 1

        print(contains_nan)

    print(empty_files)

    # for state_id in range(len(n_states)):
    #                 # params = "--input_csv /wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv --job_id {} --w_xray {} --n_state {} --init_weights rand --sa {{step3000,T300,dofA,pdb1,w1,res0}} --steps 2".format(job_id, w_xray, n_state)
    #                 # offset = str(out_dir).split("_")[-1]

    #                 # if redo:
    #                 #     bash_str = "qsub -N w{} -l h_rt={} -l mem_free=1G -l scratch=1G -t 1-1 $HOME/xray/sample_bench/scripts/sample/run_slave.sh {} {} '{}' {}".format(job_name, h_rt, job_name, job_dir, params, offset)
    #                 #     # print(bash_str)
    #                 #     # os.system(bash_str)

