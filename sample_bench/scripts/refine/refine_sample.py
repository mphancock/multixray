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
from refine import read_pdb_and_refine, read_pdb_and_refine_to_max_ff
from utility import pool_read_pdb, get_n_state_from_pdb_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_file")
    parser.add_argument("--max_ff", type=int)
    parser.add_argument("--start", type=int)
    parser.add_argument("--n_pdb", type=int)
    args = parser.parse_args()

    sample_file = Path(args.sample_file)
    analysis_dir = sample_file.parents[0]

    exp_name = analysis_dir.name
    max_ff = args.max_ff
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    ref_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/{}_ref_{}".format(exp_name, max_ff))
    ref_dir.mkdir(exist_ok=True)

    # analysis_dir = Path(Path.home(), "xray/sample_bench/data/7mhf", exp_name)
    # sample_file = Path(analysis_dir, "sample.csv")
    sample_df = pd.read_csv(sample_file, index_col=0)

    start = args.start
    end = args.start + args.n_pdb
    if end > len(sample_df):
        end = len(sample_df)
    sample_df = sample_df.iloc[start:end]
    ref_df_file = Path(analysis_dir, "ref_{}_{}.csv".format(max_ff, start))

    print(len(sample_df))
    ref_df = sample_df.copy()
    ref_df["r_free"] = np.nan
    ref_df["rmsd"] = np.nan
    ref_df["ff"] = np.nan
    ref_df["pdb"] = [Path(ref_dir, "{}.pdb".format(index)) for index in ref_df.index]

    pdb_files = [Path(pdb_file) for pdb_file in sample_df["pdb"]]
    pool_params = list()
    for index in ref_df.index:
        pdb_file = Path(sample_df.loc[index, "pdb"])
        out_pdb_file = Path(ref_dir, "{}.pdb".format(index))

        params_dict = dict()
        params_dict["pdb_file"] = pdb_file
        params_dict["out_pdb_file"] = out_pdb_file
        params_dict["max_ff"] = max_ff
        params_dict["log_file"] = None
        pool_params.append(params_dict)

    # for pool_param in pool_params:
    #     read_pdb_and_refine_to_max_ff(pool_param)

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(read_pdb_and_refine_to_max_ff, pool_params)

    for pool_result in pool_results:
        # Make sure that the pdb file existed
        print(pool_result)

    ref_df.to_csv(ref_df_file)

