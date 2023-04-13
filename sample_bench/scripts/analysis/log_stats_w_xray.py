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
import log_file_analysis


if __name__ == "__main__":
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/single_md/17_w_xray_3000_4")

    n_jobs = 102

    fields = ["xray", "xray", "xray", "rmsd", "rmsd", "rmsd", "rmsd_align", "rmsd_align", "rmsd_align"]
    stats = ["min", "mean", "std", "min", "mean", "std", "min", "mean", "std"]

    all_log_min_dict = dict()
    all_log_min_dict["job_id"] = list()
    for i in range(len(fields)):
        all_log_min_dict["{}_{}".format(fields[i], stats[i])] = list()

    for job_id in range(n_jobs):
        job_dir = Path(exp_dir, str(job_id))
        print(job_dir)
        out_dirs = list(job_dir.glob("output_*"))
        log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]
        cum_stat_dict = log_file_analysis.get_cum_stat_dict_from_logs(
            log_files=log_files,
            fields=fields,
            stats=stats,
            equil=1000
        )

        all_log_min_dict["job_id"].append(job_id)
        for i in range(len(fields)):
            entry = "{}_{}".format(fields[i], stats[i])
            all_log_min_dict[entry].append(cum_stat_dict[entry])

        print(cum_stat_dict)

    log_df = pd.DataFrame(all_log_min_dict)
    csv_file = Path(Path.home(), "xray/sample_bench/data/log_stats/17_w_xray_3000_4.csv")
    log_df.to_csv(csv_file)
