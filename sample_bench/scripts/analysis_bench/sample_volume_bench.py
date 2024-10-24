from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
import pickle
import argparse
import time
import sys

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
from params import read_job_csv


"""
Return all_group_mins which is a dataframe with 1000 rows corresponding to the 1000 groups of size n and columns corresponding to the score_field plus extra_fields. The values are the best score plus corresponding extra field values for each group.

The group contains the indices of the log_files that are to be grouped together.
"""
def get_all_group_mins(
        param_dict
):
    n = param_dict["n"]
    log_min_df = param_dict["log_min_df"]
    score_field = param_dict["score_field"]
    extra_fields = param_dict["extra_fields"]

    fields = [score_field]
    fields.extend(extra_fields)

    # Dataframe for the best score plus extra fields for each group (1000 total).
    all_group_mins_df = pd.DataFrame(index=list(range(1000)), columns=fields)

    # Get only the size n group subset.
    for i in range(1000):
        group_df = log_min_df.sample(n=n)
        min_entry = group_df.loc[group_df[score_field] == group_df[score_field].min()]
        all_group_mins_df.loc[i, score_field] = min_entry[score_field].values[0]

        for field in extra_fields:
            all_group_mins_df.loc[i, field] = min_entry[field].values[0]

    print(n)
    return n, all_group_mins_df


"""
This function performs the core sample benchmark. 1000 random groups of size n are created for n between 1 and the maximum number of logs. For each group, the best score (as evaluated by field score_field) is computed (1000 total). For all groups of size n, the mean best score as well as the variance of the best scores. The distribution of other parameters (extra_fields) for the best scoring entries from the groups can be computed as well. To avoid redudnant calculations of the best scoring structure for each log_file, a lookup table is passed as well (log_min_df).

**********
Parameters:

    log_min_df: pandas dataframe containing the minimum value for each log file. It has the format of a stat_df (columns: [field_min_0, extra_fields] and rows: [[log_file]]).

**********
Returns:

    log_stats_df: pandas dataframe with the mean and standard deviation of the best score for each group of size N. The columns are [N, (score_field)_mean, (score_field)_std, (extra_fields)_mean, (extra_fields)_std].

"""
def get_sample_volume_df(
        log_min_df,
        score_field,
        extra_fields,
        max_n
):
    log_stats_cols = list()
    log_stats_cols.append("n")

    fields = [score_field]
    fields.extend(extra_fields)
    for field in fields:
        log_stats_cols.append("{}_mean".format(field))
        log_stats_cols.append("{}_std".format(field))

    log_stats_df = pd.DataFrame(index=list(range(1,max_n+1)), columns=log_stats_cols)

    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_params = list()
    for n in range(1,max_n+1):
        params_dict = dict()
        params_dict["n"] = n
        params_dict["log_min_df"] = log_min_df
        params_dict["score_field"] = score_field
        params_dict["extra_fields"] = extra_fields
        pool_params.append(params_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    t0 = time.time()

    pool_results = list()
    # for pool_param in pool_params:
    #     pool_results.append(get_all_group_mins(pool_param))
    #     break

    pool_results = pool_obj.imap(
        get_all_group_mins,
        pool_params
    )

    for result in pool_results:
        n, all_group_mins_df = result

        for field in fields:
            field_mean = np.mean(all_group_mins_df[field])
            field_std = np.std(all_group_mins_df[field])

            log_stats_df.loc[n, "{}_mean".format(field)] = field_mean
            log_stats_df.loc[n, "{}_std".format(field)] = field_std
            log_stats_df.loc[n, "n_out"] = n

    print(time.time()-t0)

    pool_obj.close()

    return log_stats_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("--job_dir")
    # parser.add_argument("--field")
    # parser.add_argument("--bonus_fields")
    # parser.add_argument("--stat_df_file")
    # parser.add_argument("--job_csv_file")
    # job_csv_file = Path(args.job_csv_file)
    # stat_df_file = Path(args.stat_df_file)

    bonus_fields = ["ff", "rmsd"]
    stat_df_file = Path(Path.home(), "xray/sample_bench/data/7mhf/187_bench_ref_10000/sample_per_out.csv")
    job_csv_file = Path(Path.home(), "xray/sample_bench/data/params/bench.csv")

    # parser.add_argument("--out_file")
    # args = parser.parse_args()
    # print(vars(args))

    analysis_dir = stat_df_file.parents[0]
    out_file = Path(analysis_dir, "volume.csv")


    # job_dir = Path(args.job_dir)
    # out_dirs = list(job_dir.glob("output*"))

    # # First, construct the lookup table that contains for each log_file entry, the minimum value for that MD log for each of the requested fields.
    # max_frames = list()
    # log_files = list()
    # for out_dir in out_dirs:
    #     log_files.append(Path(out_dir, "log.csv"))

    # if args.bonus_fields == "":
    #     bonus_fields = list()
    # else:
    #     bonus_fields = args.bonus_fields.split(",")

    # stat_df = get_stat_df(
    #     log_file_groups=[[log_file] for log_file in log_files],
    #     main_field=args.field,
    #     main_stat="min",
    #     bonus_fields=bonus_fields,
    #     N=1,
    #     equil=1,
    #     pdb_only=True,
    #     rmsd_filter=10
    # )

    # # Drop any rows with NaN values.
    # stat_df = stat_df.dropna()
    # score_stat_df.to_csv(Path(sample_bench_dir, "stat_df_{}.csv".format(args.field)))
    volume_df = pd.DataFrame()

    for k in range(10):
        job_ids = [k*2, k*2+1, 20+k*2, 20+k*2+1]
        for job_id in job_ids:
            job_params = read_job_csv(job_csv_file=job_csv_file, job_id=job_id)
            N = job_params["N"]
            J = job_params["J"]

            stat_df = pd.read_csv(stat_df_file, index_col=0)
            stat_df = stat_df.loc[stat_df["job_id"] == job_id]
            print(len(stat_df))

            # if "rmsd_0+rmsd_1" in stat_df.columns:
            #     stat_df["rmsd_0+rmsd_1"] =  stat_df["rmsd_0+rmsd_1"]/2

            # sample_bench_bonus_fields = bonus_fields.copy()
            # if "pdb" in sample_bench_bonus_fields:
            #     sample_bench_bonus_fields.remove("pdb")

            # Next, create m random groups of n log files to compute the mean and standard deviation of the minimum field values of the log files in the group. m and n are both set to the number of total log files.
            volume_n_j_df = get_sample_volume_df(
                log_min_df=stat_df,
                score_field="field",
                extra_fields=["r_free_0"],
                max_n=len(stat_df)
            )
            volume_n_j_df["k"] = k
            volume_n_j_df["job_id"] = job_id
            volume_n_j_df["N"] = N
            volume_n_j_df["J"] = J
            volume_df = pd.concat([volume_df, volume_n_j_df])

    print(volume_df.head())

    volume_df.reset_index(drop=True, inplace=True)
    volume_df.to_csv(Path(out_file))