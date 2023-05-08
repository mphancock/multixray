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
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/46_w_xray")
    log_file = Path(Path.home(), "xray/sample_bench/data/3ca7/46_w_xray/r_free.csv")
    n_jobs = 80

    log_df = pd.DataFrame(index=list(range(n_jobs)), columns=["avg_min_r_free", "min_r_free"])

    for job_id in range(n_jobs+1):
        job_dir = Path(exp_dir, str(job_id))
        print(job_dir)

        out_dirs = list(job_dir.glob("output_*"))

        # Get the average min r free and min r free from all output logs for a given w_xray.
        log_file_groups = [[Path(out_dir, "log.csv")] for out_dir in out_dirs]
        stat_df = get_stat_df.get_stat_df(
            log_file_groups=log_file_groups,
            fields=["r_free"],
            stats=["min"],
            equil=1000
        )
        stat_df.dropna(inplace=True)

        min_r_free = stat_df["r_free_min"].min()
        avg_min_r_free = stat_df["r_free_min"].mean()

        log_df.loc[job_id, "avg_min_r_free"] = avg_min_r_free
        log_df.loc[job_id, "min_r_free"] = min_r_free

    log_df.to_csv(log_file)