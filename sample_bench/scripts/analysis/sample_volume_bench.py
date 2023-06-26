from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
import pickle
import argparse
import time
import sys

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df
import sample_bench
sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--job_dir")
    parser.add_argument("--sample_bench_dir")
    parser.add_argument("--ref_pdb_file")
    parser.add_argument("--field")
    parser.add_argument("--bonus_fields")
    parser.add_argument("--max_n", type=int)
    args = parser.parse_args()
    print(args.job_dir)
    print(args.sample_bench_dir)
    print(args.ref_pdb_file)
    print(args.field)
    print(args.bonus_fields)
    print(args.max_n)

    job_dir = Path(args.job_dir)
    out_dirs = list(job_dir.glob("output*"))
    out_dirs = out_dirs

    sample_bench_dir = Path(args.sample_bench_dir)
    sample_bench_dir.mkdir(parents=True, exist_ok=True)

    # First, construct the lookup table that contains for each log_file entry, the minimum value for that MD log for each of the requested fields.
    max_frames = list()
    log_files = list()
    for out_dir in out_dirs:
        log_files.append(Path(out_dir, "log.csv"))

    bonus_fields = args.bonus_fields.split(",")
    score_field = args.field
    score_stat_df = get_stat_df.get_stat_df(
        log_file_groups=[[log_file] for log_file in log_files],
        main_field=score_field,
        main_stat="min",
        bonus_fields=bonus_fields,
        N=1,
        equil=100,
        pdb_only=True
    )
    score_stat_df.to_csv(Path(sample_bench_dir, "stat_df_{}.csv".format(score_field)))
    # score_stat_df = pd.read_csv(Path(sample_bench_dir, "stat_df_{}.csv".format(score_field)), index_col=0)

    if args.max_n > len(log_files):
        max_n = len(log_files)
    else:
        max_n = args.max_n

    sample_bench_bonus_fields = bonus_fields.copy()
    if "pdb" in sample_bench_bonus_fields:
        sample_bench_bonus_fields.remove("pdb")

    # Next, create m random groups of n log files to compute the mean and standard deviation of the minimum field values of the log files in the group. m and n are both set to the number of total log files.
    sample_volume_df = sample_bench.get_sample_volume_df(
        log_files=log_files,
        log_min_df=score_stat_df,
        score_field=score_field,
        extra_fields=sample_bench_bonus_fields,
        max_n=max_n
    )
    sample_volume_df.to_csv(Path(sample_bench_dir, "sample_bench_{}.csv".format(score_field)))