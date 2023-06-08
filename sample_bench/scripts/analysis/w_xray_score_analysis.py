from pathlib import Path
import sys
import multiprocessing
import pandas as pd
import os
import numpy as np

import IMP
import IMP.atom
import IMP.core

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df


if __name__ == "__main__":
    job_name = "40_wxray_synth"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", job_name)
    analysis_dir = Path(Path.home(), "xray/sample_bench/data/7mhf", job_name)
    analysis_dir.mkdir(exist_ok=True)
    log_file = Path(analysis_dir, "r_free.csv".format(job_name))
    n_jobs = 70

    log_df = pd.DataFrame(index=list(range(n_jobs)), columns=["avg_min_r_free", "min_r_free"])

    for job_id in range(n_jobs+1):
        job_dir = Path(exp_dir, str(job_id))
        print(job_dir)
        out_dirs = list(job_dir.glob("output_*"))
        log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]

        # Get the average min R free and min R free from all output logs for a given w_xray.
        # R_free_0 is the first R_free reported.
        field = "r_free_0"
        stat_df = get_stat_df.get_stat_df(
            log_file_groups=[[log_file] for log_file in log_files],
            fields=[field],
            stats=["min"],
            N=1,
            offset=1,
            equil=1000
        )
        stat_df.dropna(inplace=True)
        min_r_free = stat_df["{}_min_0".format(field)].min()
        avg_min_r_free = stat_df["{}_min_0".format(field)].mean()

        # Get the mean R free across all logs for a given w_xray.
        mean_stat_df = get_stat_df.get_stat_df(
            log_file_groups=[log_files],
            fields=[field],
            stats=["mean"],
            N=1,
            offset=1,
            equil=1000
        )
        mean_stat_df.dropna(inplace=True)
        avg_r_free = mean_stat_df.iloc[0]["{}_mean".format(field)]

        log_df.loc[job_id, "avg_r_free"] = avg_r_free
        log_df.loc[job_id, "avg_min_r_free"] = avg_min_r_free
        log_df.loc[job_id, "min_r_free"] = min_r_free

    log_df.to_csv(log_file)