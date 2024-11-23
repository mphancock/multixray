from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
import pickle
import argparse
import time
import sys

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
            log_stats_df.loc[n, "{}_50".format(field)] = np.percentile(all_group_mins_df[field], 50)
            log_stats_df.loc[n, "{}_25".format(field)] = np.percentile(all_group_mins_df[field], 25)
            log_stats_df.loc[n, "{}_75".format(field)] = np.percentile(all_group_mins_df[field], 75)
            log_stats_df.loc[n, "n"] = n

    print(time.time()-t0)

    pool_obj.close()

    return log_stats_df


if __name__ == "__main__":
    # field = "xray_native_0+xray_native_1"
    # bonus_fields = ["ff", "rmsd"]
    field = "rmsd"
    bonus_fields = []

    stat_df_file = Path(Path.home(), "xray/sample_bench/data/analysis/267_full_ref/all_outs_0.csv")
    job_csv_file = Path(Path.home(), "xray/sample_bench/data/params/267.csv")
    job_id = 0

    analysis_dir = stat_df_file.parents[0]
    out_file = Path(analysis_dir, "volume_{}.csv".format(field))

    volume_df = pd.DataFrame()

    job_params = read_job_csv(job_csv_file=job_csv_file, job_id=job_id)
    N = job_params["N"]
    J = job_params["J"]

    stat_df = pd.read_csv(stat_df_file, index_col=0)
    print(len(stat_df))

    # Next, create m random groups of n log files to compute the mean and standard deviation of the minimum field values of the log files in the group. m and n are both set to the number of total log files.
    volume_df = get_sample_volume_df(
        log_min_df=stat_df,
        score_field=field,
        extra_fields=bonus_fields,
        max_n=len(stat_df)
    )

    print(volume_df.head())

    volume_df.reset_index(drop=True, inplace=True)
    volume_df.to_csv(Path(out_file))