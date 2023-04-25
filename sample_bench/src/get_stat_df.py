from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing


"""
get_stat_from_log_files takes a list of log_files and returns a statistic (eg, mean) of a given field post equilibration.

log_files: list, list of log files.
equil: int, number of lines to skip in the log file.
field: str, name of the field to compute the statistic on.
stat: str, name of the statistic to compute.
"""
def pool_get_stat_from_log_files(
        params_dict
):
    log_files = params_dict["log_files"]
    equil = params_dict["equil"]
    field = params_dict["field"]
    stat = params_dict["stat"]

    # print(log_files, equil, field, stat)
    log_files, field, stat, log_stat, stat_id = get_stat_from_log_files(
        log_files=log_files,
        equil=equil,
        field=field,
        stat=stat
    )

    return log_files, field, stat, log_stat, stat_id


def get_stat_from_log_files(
        log_files,
        equil,
        field,
        stat
):
    # Merge the log files into a single dataframe.
    log_dfs = list()
    for log_file in log_files:
        try:
            # print("trying to read log file: {}".format(log_file))
            log_df = pd.read_csv(log_file)
        except pd.errors.EmptyDataError:
            # print("catching empty data error")
            continue
        log_dfs.append(log_df[equil:])

    merege_log_df = pd.concat(log_dfs)

    field, stat, log_stat, stat_id = get_stat_from_log_df(
        log_df=merege_log_df,
        field=field,
        stat=stat
    )

    return log_files, field, stat, log_stat, stat_id


"""
get_stat_from_log_df takes a log_df and returns a statistic (eg, mean) of a given field post equilibration.

log_df: pandas dataframe, log file as a dataframe.
field: str, name of the field to compute the statistic on.
stat: str, name of the statistic to compute.
"""
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


"""
get_stat_df function takes in a set of output directories (presumably the set of output directories for a job) and returns a dataframe containing the stat for each field for each log file for all stat, field pair. Having this dataframe is important for sample volume benchmarks.

log_file_groups: the set of groups where each group contains one or more log files (the rows of the stat_df). Each stat, field pair is computed for each group of log files.
fields: the fields (eg, rmsd) from each log file to extract a statistic for.
stats: the statistic (eg, mean) to compute for each field (field[i]_statistic[i] is a column).
equil: the number of frames to discard before computing a statistic.
include_id: whether to include an additional column for each requested min/max which records the relevant file.
"""
def get_stat_df(
        log_file_groups,
        fields,
        stats,
        equil,
        include_id
):
    for stat in stats:
        if stat not in ["min", "max", "mean", "std", "var"]:
            raise RuntimeError("Invalid stat selection: {}".format(stat))

    n_fields = len(fields)

    # Construct the columns for the stat_df.
    columns = list()
    for i in range(n_fields):
        column = "{}_{}".format(fields[i], stats[i])
        columns.append(column)
        if stats[i] in ["min", "max"] and include_id:
            columns.append("{}_id".format(column))

    indices = [str(group) for group in log_file_groups]
    log_stat_df = pd.DataFrame(index=indices, columns=columns)

    # Iterate through each group of log files to distribute the calculation of a each field, stat pair.
    pool_params = list()
    for log_file_group in log_file_groups:
        # For each field, stat pair, create a pool_param.
        for i in range(n_fields):
            pool_param = dict()
            pool_param["log_files"] = log_file_group
            pool_param["equil"] = equil
            pool_param["field"] = fields[i]
            pool_param["stat"] = stats[i]
            pool_params.append(pool_param)

    # print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_get_stat_from_log_files,
        pool_params
    )

    # Iterate through the pool results and fill in the stat_df.
    # print(log_stat_df.head())
    for pool_result in pool_results:
        log_files, field, stat, log_stat, stat_id = pool_result
        column = "{}_{}".format(field, stat)
        entry = str(log_files)
        log_stat_df.loc[entry, column] = log_stat

        if stat_id and include_id:
            log_stat_df.loc[entry, log_stat_df.columns.get_loc("{}_id".format(column))] = stat_id

    pool_obj.close()
    return log_stat_df


# if __name__ == "__main__":
#     job_dir_1 = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/46_w_xray/0")
#     job_dir_2 = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/46_w_xray/1")

#     for log_file in job_dir.glob("output_*/log.csv"):
#         log_file_groups = list()
#         log_file_groups.append([log_file])

#         log_stat_df = get_stat_df(
#             log_file_groups=log_file_groups,
#             fields=["xray"],
#             stats=["mean"],
#             equil=1000,
#             include_id=True
#         )

#         print(log_stat_df.head())

#     log_file_groups = [tuple(job_dir_1.glob("output_*/log.csv")), tuple(job_dir_2.glob("output_*/log.csv"))]
#     log_stat_df = get_stat_df(
#         log_file_groups=log_file_groups,
#         fields=["xray"],
#         stats=["mean"],
#         equil=1000,
#         include_id=True
#     )

#     print(log_stat_df.head())