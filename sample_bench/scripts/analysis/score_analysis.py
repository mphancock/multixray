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
    job_name = "66_native_1x"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7", job_name)
    analysis_dir = Path(Path.home(), "xray/sample_bench/data/3ca7", job_name)
    analysis_dir.mkdir(exist_ok=True)
    log_file = Path(analysis_dir, "score_analysis.csv".format(job_name))
    n_jobs = 40
    equil = 1

    log_df = pd.DataFrame(index=list(range(n_jobs)))

    # for job_id in range(n_jobs+1):
    for job_id in range(n_jobs):
        job_dir = Path(exp_dir, str(job_id))
        print(job_dir)
        out_dirs = list(job_dir.glob("output_*"))
        log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]

        # Get the average min R free and min R free from all output logs for a given w_xray.
        # R_free_0 is the first R_free reported.
        field = "xray_0"
        bonus_fields = ["r_free_0", "rmsd_ord", "rmsd_avg", "pdb"]
        stat_df = get_stat_df.get_stat_df(
            log_file_groups=[[log_file] for log_file in log_files],
            main_field=field,
            main_stat="min",
            bonus_fields=bonus_fields,
            N=1,
            equil=1,
            pdb_only=True,
            test=False
        )
        print(len(stat_df))
        stat_df.to_csv(Path(Path.home(), "xray/tmp/stat_df.csv"))

        stat_df.dropna(inplace=True)

        print(stat_df.columns)

        log_df.loc[job_id, "min_xray"] = stat_df["{}_min_0".format(field)].min()
        log_df.loc[job_id, "avg_min_xray"] = stat_df["{}_min_0".format(field)].mean()

        for bonus_field in bonus_fields:
            if bonus_field == "pdb":
                continue

            log_df.loc[job_id, "min_xray_{}".format(bonus_field)] = stat_df["{}_min_0_{}".format(field, bonus_field)].min()
            log_df.loc[job_id, "avg_min_xray_{}".format(bonus_field)] = stat_df["{}_min_0_{}".format(field, bonus_field)].mean()

        print(len(stat_df))

        stat_df["xray_0_min_0"] = pd.to_numeric(stat_df["xray_0_min_0"])
        print(stat_df["xray_0_min_0"].min())

        print(stat_df["xray_0_min_0"].idxmin())
        log_df.loc[job_id, "pdb"] = stat_df.loc[stat_df["{}_min_0".format(field)].idxmin(), "xray_0_min_0_pdb"]

    log_df.to_csv(log_file)