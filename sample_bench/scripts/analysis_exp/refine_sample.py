from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import argparse

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import read_pdb_and_refine
from utility import pool_read_pdb, get_n_state_from_pdb_file
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis")))
from score_analysis_exp import field_id_mapper, map_std_field_to_field



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp_name")
    parser.add_argument("--n_step", type=int)
    args = parser.parse_args()

    exp_name = args.exp_name
    n_step = args.n_step
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    ref_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/{}_ref_{}".format(exp_name, n_step))
    ref_dir.mkdir(exist_ok=True)

    sample_df = pd.read_csv(Path(Path.home(), "xray/sample_bench/data/7mhf/{}/sample.csv".format(exp_name)), index_col=0)

    print(len(sample_df))

    # print(len(stat_df))
    # stat_df.to_csv(Path(Path.home(), "xray/dev/36_refine_sample/data/stat_df.csv"))

    analysis_dir = Path(Path.home(), "xray/sample_bench/data/7mhf", exp_name)

    ref_df = sample_df.copy()
    ref_df["pdb"] = [Path(ref_dir, "{}.pdb".format(i)) for i in range(len(ref_df))]
    ref_df.to_csv(Path(analysis_dir, "ref_{}.csv".format(n_step)))

    pdb_files = [Path(pdb_file) for pdb_file in sample_df["pdb"]]
    pool_params = list()
    pdb_id = 0
    for pdb_file in pdb_files:
        params_dict = dict()
        params_dict["pdb_file"] = pdb_file
        params_dict["out_pdb_file"] = Path(ref_dir, "{}.pdb".format(pdb_id))
        params_dict["n_step"] = n_step
        params_dict["log_file"] = None
        pool_params.append(params_dict)

        pdb_id = pdb_id + 1

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(read_pdb_and_refine, pool_params)

    for pool_result in pool_results:
        # Make sure that the pdb file existed
        print(pool_result)

