from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
import pickle
import argparse
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


# The function takes a log and returns a statistic (eg, mean) of a given field post equilibration.
def pool_get_stat_from_log_df(
        params_dict
):
    log_file = params_dict["log_file"]
    field, stat, log_stat, stat_id = get_stat_from_log_df(
        log_df=params_dict["log_df"],
        field=params_dict["field"],
        stat=params_dict["stat"]
    )
    return log_file, field, stat, log_stat, stat_id


def get_stat_from_log_df(
        log_df,
        field,
        stat,
):
    if field == "tot":
        log_df["tot"] = 30000*log_df["xray"] + log_df["ff"]

    stat_id = None
    if stat == "min":
        log_stat = log_df[field].min()
        stat_id = log_df[field].idxmin()
    elif stat == "max":
        log_stat = log_df[field].max()
        stat_id = log_df[field].idxmax()
    elif stat == "mean":
        log_stat = log_df[field].mean()
    elif stat == "std":
        log_stat = log_df[field].std()
    elif stat == "var":
        log_stat = log_df[field].var()
    else:
        raise RuntimeError("Invalid stat selection: {}".format(stat))

    return field, stat, log_stat, stat_id



# This function is similar to get_stat_df_from_logs except for each field and stat it returns the field stat for all log files put together (cumulative) as a dictionary. This is useful for the w_xray analysis.
def get_cum_stat_dict_from_logs(
        log_files,
        fields,
        stats,
        equil
):
    pool_params = list()
    n_fields = len(fields)

    log_dfs = list()
    for log_file in log_files:
        try:
            log_df = pd.read_csv(log_file, index_col=0)
        except pd.errors.EmptyDataError:
            continue

        if len(log_df) < equil:
            continue

        log_df_equil = log_df.iloc[equil:]
        log_dfs.append(log_df_equil)

    merge_log_file_df = pd.concat(log_dfs)

    for i in range(n_fields):
        pool_param = dict()
        pool_param["log_file"] = None
        pool_param["log_df"] = merge_log_file_df
        pool_param["field"] = fields[i]
        pool_param["stat"] = stats[i]
        pool_params.append(pool_param)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_get_stat_from_log_df,
        pool_params
    )

    cum_stat_dict = dict()
    for pool_result in pool_results:
        if not pool_result:
            continue
        else:
            log_file, field, stat, log_stat, stat_pdb_file = pool_result

        cum_stat_dict["{}_{}".format(field, stat)] = log_stat

    return cum_stat_dict


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

    return log_stats_df
