from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import os

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from stat_df import get_stat_df
# from refine import refine, pool_refine
from utility import pool_read_pdb
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis_exp")))


def check_job_dir(
    pool_params
):
    job_dir = pool_params["job_dir"]
    job_num = job_dir.name
    n_out_dir = pool_params["n_out_dir"]
    job_csv_file = pool_params["job_csv_file"]

    job_name = job_dir.name
    # job_name = "n{}_j{}".format(state_id, job_id)
    # job_dir = Path(exp_dir, job_name)
    out_dirs = [Path(job_dir, "output_{}".format(i)) for i in range(n_out_dir)]

    n_valid = 0
    avg_n_pdb_files = 0
    run_times = list()
    for out_dir in out_dirs:
        out_dir_name = out_dir.name
        valid = True
        if not out_dir.exists():
            print("{} {} does not exist".format(job_num, out_dir_name))
            valid = False
            n_pdb_files = 0
            # continue
        else:
            if not Path(out_dir, "log.csv").exists():
                print("{} {} does not contain log.csv".format(job_num, out_dir_name))
                valid = False
                n_pdb_files = 0
            else:
                pdb_files = out_dir.glob("pdbs/*")
                n_pdb_files = len(list(pdb_files))

                if n_pdb_files < 250:
                    print("{} {} does not contain enough ({})".format(job_num, out_dir_name, n_pdb_files))
                    valid = False

        avg_n_pdb_files += n_pdb_files

    if len(out_dirs) > 0:
        avg_n_pdb_files /= len(out_dirs)

    if n_valid > 0:
        run_time_avg = np.mean(run_times)
        run_time_std = np.std(run_times)
    else:
        run_time_avg = None
        run_time_std = None

    return job_name, n_valid, avg_n_pdb_files, run_time_avg, run_time_std


def get_run_time(
    log_file
):
    log_df = pd.read_csv(log_file)
    get_run_time = log_df.loc[len(log_df)-1, "time"]

    return get_run_time


if __name__ == "__main__":
    exp_name = "280_exp_all_2"
    job_csv_file = "/wynton/home/sali/mhancock/xray/sample_bench/data/params/280.csv"
    n_out_dir = 1000

    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out", exp_name)
    analysis_dir = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/analysis", exp_name)
    # analysis_dir.mkdir(exist_ok=True)

    job_stats_df = pd.DataFrame()

    job_dirs = exp_dir.glob("*")
    params = list()
    for job_dir in job_dirs:
        exp_dir = job_dir.parents[0]

        params_dict = dict()
        params_dict["job_dir"] = job_dir
        params_dict["n_out_dir"] = n_out_dir
        params_dict["job_csv_file"] = job_csv_file
        params.append(params_dict)

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(check_job_dir, params)

    for pool_result in pool_results:
        job_name, n_valid, avg_n_pdb_files, run_time_avg, run_time_std = pool_result

        df_id = len(job_stats_df)
        job_stats_df.loc[df_id, "job_name"] = job_name
        # job_stats_df.loc[df_id, "n_out_dir"] = len(out_dirs)
        job_stats_df.loc[df_id, "n_valid"] = n_valid
        job_stats_df.loc[df_id, "avg_n_pdb_files"] = avg_n_pdb_files

        print(job_name, n_valid, avg_n_pdb_files, run_time_avg, run_time_std)

