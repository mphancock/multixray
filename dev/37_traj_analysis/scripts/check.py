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
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis_exp")))


def get_run_time(
    log_file
):
    log_df = pd.read_csv(log_file)
    get_run_time = log_df.loc[len(log_df)-1, "time"]

    return get_run_time


if __name__ == "__main__":
    exp_name = "186_w_xray_bench"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    analysis_dir = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/7mhf", exp_name)
    analysis_dir.mkdir(exist_ok=True)

    n_out_dir = 5
    h_rt = "8:00:00"

    n_states = [1, 2, 4, 8, 16]

    job_stats_df = pd.DataFrame()

    job_dirs = exp_dir.glob("*")
    for job_dir in job_dirs:
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
                print("{} does not exist".format(out_dir_name))
                valid = False
            else:
                pdb_files = out_dir.glob("pdbs/*")
                n_pdb_files = len(list(pdb_files))

                if not Path(out_dir, "log.csv").exists():
                    print("{} does not contain log.csv".format(out_dir_name))
                    valid = False
                elif n_pdb_files < 302:
                    print("{} does not contain enough ({})".format(out_dir_name, n_pdb_files))
                    valid = False

            if valid:
                n_valid += 1
                run_time = get_run_time(Path(out_dir, "log.csv"))
                run_times.append(run_time)
            else:
                offset = str(out_dir).split("_")[-1]

                params = "--params_file {}".format(Path(out_dir, "params.csv"))

                bash_str = "qsub -N b{} -l h_rt={} -l mem_free=1G -l scratch=1G -t 1-1 /wynton/home/sali/mhancock/xray/sample_bench/scripts/sample/run_slave.sh {} {} '{}'".format(job_name, h_rt, job_dir, offset, params)

            #     print(bash_str)
            # #     print(bash_str)
            #     os.system(bash_str)


            #     # params = "--input_csv /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/csvs/7mhf_20.csv --job_id {} --w_xray XX --n_state {} --init_weights rand --sa {{step3000,T300,dofA,pdb1,w1,res1.5}} --steps 2".format(job_id, n_state)

            #     # bash_str = "qsub -N w{} -l h_rt={} -l mem_free=1G -l scratch=1G -t 1-1 $HOME/xray/sample_bench/scripts/sample/run_slave.sh {} {} '{}' {}".format(job_name, h_rt, job_name, job_dir, params, offset)

            #     orig_exp_name = exp_name.split("_ref_")[0]
            #     orig_job_dir = Path(exp_dir.parents[0], orig_exp_name, job_name)


            avg_n_pdb_files += n_pdb_files

        if len(out_dirs) > 0:
            avg_n_pdb_files /= len(out_dirs)

        if n_valid > 0:
            run_time_avg = np.mean(run_times)
            run_time_std = np.std(run_times)
        else:
            run_time_avg = None
            run_time_std = None

        df_id = len(job_stats_df)
        job_stats_df.loc[df_id, "job_name"] = job_name
        job_stats_df.loc[df_id, "n_out_dir"] = len(out_dirs)
        job_stats_df.loc[df_id, "n_valid"] = n_valid
        job_stats_df.loc[df_id, "avg_n_pdb_files"] = avg_n_pdb_files

        print(job_dir.name, len(out_dirs), n_valid, avg_n_pdb_files, run_time_avg, run_time_std)

