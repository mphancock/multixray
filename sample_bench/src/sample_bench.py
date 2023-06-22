from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
import pickle
import time
import sys


"""
Use the lookup table to get the best score (as evaluated by score_field) for each log file group in addition to the corresponding extra fields for that best scoring frame.
"""
def get_field_min_for_group(
        log_file_group,
        log_min_df,
        score_field,
        extra_fields
):
    group_mins = dict()
    fields = [score_field]
    fields.extend(extra_fields)
    for field in fields:
        group_mins[field] = np.inf

    for log_file in log_file_group:
        log_file_index = [str(log_file)]
        score_field_column = "{}_min_0".format(score_field)
        log_min = float(log_min_df.loc[log_file_index, score_field_column])

        if log_min < group_mins[score_field]:
            group_mins[score_field] = log_min

            for field in extra_fields:
                extra_field_column = "{}_{}".format(score_field_column, field)
                group_mins[field] = float(log_min_df.loc[log_file_index, extra_field_column])

    return group_mins


"""
Return all_group_mins which is a dataframe with 1000 rows corresponding to the 1000 groups of size n and columns corresponding to the score_field plus extra_fields. The values are the best score plus corresponding extra field values for each group.

The group contains the indices of the log_files that are to be grouped together.
"""
def get_all_group_mins(
        param_dict
):
    n = param_dict["n"]
    # groups_file = param_dict["groups_file"]
    log_files = param_dict["log_files"]
    log_min_df = param_dict["log_min_df"]
    score_field = param_dict["score_field"]
    extra_fields = param_dict["extra_fields"]

    groups_file = Path(Path.home(), "xray/sample_bench/data/etc/1000.pl")

    with open(groups_file, 'rb') as f:
        all_groups_1000 = pickle.load(f)

    fields = [score_field]
    fields.extend(extra_fields)

    # Dataframe for the best score plus extra fields for each group (1000 total).
    all_group_mins_df = pd.DataFrame(index=list(range(1000)), columns=fields)

    # Get only the size n group subset.
    all_groups_n = list()
    for group_1000 in all_groups_1000:
        group_n = group_1000[:n]
        all_groups_n.append(group_n)

    # Group mins contains the minimum value for each group.
    for i in range(len(all_groups_n)):
        group_ids = all_groups_n[i]

        # Need to handle the case where the ids exceed the size of log_files.
        log_file_group = list()
        for id in group_ids:
            if id < len(log_files):
                log_file_group.append(log_files[id])
            else:
                continue

        # Dictionary containing best score for a group plus corresponding extra field values.
        group_mins = get_field_min_for_group(
            log_file_group=log_file_group,
            log_min_df=log_min_df,
            score_field=score_field,
            extra_fields=extra_fields
        )

        for field in fields:
            all_group_mins_df.loc[i, field] = group_mins[field]

    print(n)
    return n, all_group_mins_df


"""
This function performs the core sample benchmark. 1000 random groups of size n are created for n between 1 and the maximum number of logs. For each group, the best score (as evaluated by field score_field) is computed (1000 total). For all groups of size n, the mean best score as well as the variance of the best scores. The distribution of other parameters (extra_fields) for the best scoring entries from the groups can be computed as well. To avoid redudnant calculations of the best scoring structure for each log_file, a lookup table is passed as well (log_min_df).

Parameters:
    log_files: list of Path objects to the log files.
    log_min_df: pandas dataframe containing the minimum value for each log file. It has the format of a stat_df (columns: [field_min_0, extra_fields] and rows: [[log_file]])

"""
def get_sample_volume_df(
        log_files,
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
        params_dict["log_files"] = log_files
        params_dict["log_min_df"] = log_min_df
        params_dict["score_field"] = score_field
        params_dict["extra_fields"] = extra_fields
        pool_params.append(params_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    t0 = time.time()
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

    print(time.time()-t0)

    pool_obj.close()

    return log_stats_df