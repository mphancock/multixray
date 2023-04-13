from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
import pickle
import argparse
import time
import sys

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import log_file_analysis


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--job_name")
    parser.add_argument("--job_id", type=int)
    parser.add_argument("--max_n", type=int)
    args = parser.parse_args()
    print(args.job_name)
    print(args.job_id)
    print(args.max_n)

    job_name = args.job_name
    exp_id = job_name.split("_")[0]
    job_id = args.job_id

    xray_dir = Path("/wynton/group/sali/mhancock/xray")
    out_dirs = list(Path(xray_dir, "sample_bench/out/3ca7/single_md/{}/{}".format(job_name, job_id)).glob("output*"))
    fields = ["tot", "rmsd"]

    # First, construct the lookup table that contains for each log_file entry, the minimum value for that MD log for each of the requested fields.
    log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]
    stat_df = log_file_analysis.get_stat_df_from_logs(
        log_files=log_files,
        fields=["tot", "rmsd"],
        stats=["min", "min"],
        equil=0,
        include_id=True
    )
    print(stat_df.head())

    # Create an additional column that is the rmsd of the minimum scoring structure.
    min_score_rmsds = list()
    min_score_rmsd_aligns = list()
    for i in range(len(stat_df)):
        log_file = stat_df.loc[i, "log_file"]
        min_score_id = stat_df.loc[i, "tot_min_id"]

        log_df = pd.read_csv(log_file, index_col=0)
        min_score_rmsd = log_df.loc[min_score_id, "rmsd"]
        min_score_rmsd_align = log_df.loc[min_score_id, "rmsd_align"]
        min_score_rmsds.append(min_score_rmsd)
        min_score_rmsd_aligns.append(min_score_rmsd_align)

    stat_df["tot_min_rmsd"] = min_score_rmsds
    stat_df["tot_min_rmsd_align"] = min_score_rmsd_aligns

    stat_df.to_csv(Path(Path.home(), "xray/sample_bench/data/stat_df/{}.csv".format(exp_id)))

    # stat_df = pd.read_csv(Path(Path.home(), "xray/sample_bench/data/stat_df/{}.csv".format(exp_id)), index_col=0)
    # for i in range(len(stat_df)):
    #     stat_df.loc[i, "log_file"] = Path(stat_df.loc[i, "log_file"])
    #
    # if args.max_n > len(log_files):
    #     max_n = len(log_files)
    # else:
    #     max_n = args.max_n
    #
    # # Next, create m random groups of n log files to compute the mean and standard deviation of the minimum field values of the log files in the group. m and n are both set to the number of total log files.
    # sample_volume_df = log_file_analysis.get_sample_volume_df(
    #     log_files=log_files,
    #     log_min_df=stat_df,
    #     fields=["tot_min", "rmsd_min", "tot_min_rmsd"],
    #     max_n=max_n
    # )
    # sample_volume_df.to_csv(Path(Path.home(), "xray/sample_bench/data/volume_bench/{}.csv".format(exp_id)))
    #

