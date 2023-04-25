from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
import pickle
import time


def get_field_min_from_log_files(
        log_files,
        log_min_df,
        field
):
    group_min_score = 10**10

    for log_file in log_files:
        log_min = log_min_df[log_min_df["log_file"] == log_file][field].min()
        if log_min < group_min_score:
            group_min_score = log_min

    return group_min_score


def get_all_group_mins(
        param_dict
):
    n = param_dict["n"]
    groups_file = param_dict["groups_file"]
    log_files = param_dict["log_files"]
    log_min_df = param_dict["log_min_df"]
    field = param_dict["field"]

    with open(groups_file, 'rb') as f:
        all_groups_size_n = pickle.load(f)

    all_group_mins = list()
    for group_ids in all_groups_size_n:
        # Need to handle the case where the ids exceed the size of log_files.
        group = list()
        # group = [log_files[id] for id in ids_group]
        for id in group_ids:
            if id < len(log_files):
                group.append(log_files[id])
            else:
                continue

        group_min = get_field_min_from_log_files(
            log_files=group,
            log_min_df=log_min_df,
            field=field
        )

        # Must remove large group minimums that indicate numerical instability within the simulation.
        if group_min < 10**6:
            all_group_mins.append(group_min)

    return n, field, all_group_mins


def get_sample_volume_df(
        log_files,
        log_min_df,
        fields,
        max_n
):
    groups_dir = Path(Path.home(), "xray/sample_bench/data/groups/m_n_1000")

    log_stats_cols = list()
    log_stats_cols.append("n")
    for field in fields:
        log_stats_cols.append("{}_mean".format(field))
        log_stats_cols.append("{}_std".format(field))

    log_stats_df = pd.DataFrame(0, index=list(range(max_n)), columns=log_stats_cols)
    log_stats_df["n"] = list(range(1,max_n+1))

    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_params = list()
    for n in range(1,max_n+1):
        # if n < max_n:
        #     continue

        groups_file = Path(groups_dir, "{}.pl".format(n))
        for field in fields:
            params_dict = dict()
            params_dict["n"] = n
            params_dict["groups_file"] = groups_file
            params_dict["log_files"] = log_files
            params_dict["log_min_df"] = log_min_df
            params_dict["field"] = field
            pool_params.append(params_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    t0 = time.time()
    pool_results = pool_obj.imap(
        get_all_group_mins,
        pool_params
    )

    for result in pool_results:
        n, field, all_group_mins = result
        groups_size_n_mean = np.mean(all_group_mins)
        groups_size_n_std = np.std(all_group_mins)
        print(n, field, groups_size_n_mean, groups_size_n_std)

        log_stats_df.iloc[n-1, log_stats_df.columns.get_loc("{}_mean".format(field))] = groups_size_n_mean
        log_stats_df.iloc[n-1, log_stats_df.columns.get_loc("{}_std".format(field))] = groups_size_n_std
    print(time.time()-t0)

    pool_obj.close()

    return log_stats_df

