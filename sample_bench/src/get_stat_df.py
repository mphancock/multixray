from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing


"""
get_stat_from_log_df takes a log_df and returns a statistic (eg, mean) of a given field post equilibration.

params
log_df: pandas dataframe, log file as a dataframe.
field: str, name of the field to compute the statistic on.
stat: str, name of the statistic to compute.
N: int, the number of values to return for a stat. Only applicable to min and max.

returns
field_stat_df: pandas dataframe, dataframe containing the stat value(s) plus entry information. The dataframe will only be greater than a single row for min and max with N>1.
"""
def get_stat_from_log_df(
        log_df,
        field,
        stat,
        N
):
    field_stat_df = pd.DataFrame(columns=["field", "stat", "value", "log", "id"])

    # 3 error cases to eturn nan entries: if the field is not in the log_df, if the log_df is empty, or if the stat is not invalid.
    if field not in log_df.columns:
        field_stat_df.loc[0] = [field, stat, np.nan, np.nan, np.nan]
    elif len(log_df) == 0:
        field_stat_df.loc[0] = [field, stat, np.nan, np.nan, np.nan]
    else:
        if stat in ["min", "max"]:
            if stat == "min":
                n_stat_df = log_df.nsmallest(N, field)
            else:
                n_stat_df = log_df.nlargest(N, field)

            for i in range(len(n_stat_df)):
                field_stat_df.loc[i] = [field, stat, n_stat_df[field].iloc[i], n_stat_df["log"].iloc[i], n_stat_df["id"].iloc[i]]
        elif stat == "mean":
            field_stat_df.loc[0] = [field, stat, log_df[field].mean(), np.nan, np.nan]
        elif stat == "std":
            field_stat_df.loc[0] = [field, stat, log_df[field].std(), np.nan, np.nan]
        elif stat == "var":
            field_stat_df.loc[0] = [field, stat, log_df[field].var(), np.nan, np.nan]
        else:
            field_stat_df.loc[0] = [field, stat, np.nan, np.nan, np.nan]

    return field_stat_df


"""
get_stat_from_log_files takes a list of log_files and returns a statistic (eg, mean) of a given field post equilibration.

params
log_files: list, list of log files.
equil: int, number of lines to skip in the log file.
field: str, name of the field to compute the statistic on.
stat: str, name of the statistic to compute.

returns
field_stat_df: pandas dataframe, dataframe containing the stat value(s) plus entry information. The dataframe will only be greater than a single row for min and max with N>1.
"""
def pool_get_stat_from_log_files(
        params_dict
):
    log_files = params_dict["log_files"]
    equil = params_dict["equil"]
    field = params_dict["field"]
    stat = params_dict["stat"]
    N = params_dict["N"]
    offset = params_dict["offset"]

    # print(log_files, equil, field, stat)
    log_files, field_stat_df = get_stat_from_log_files(
        log_files=log_files,
        equil=equil,
        field=field,
        stat=stat,
        N=N,
        offset=offset
    )

    return log_files, field_stat_df


def get_stat_from_log_files(
        log_files,
        field,
        stat,
        equil,
        N,
        offset
):
    # Merge the log files into a single dataframe.
    log_dfs = list()
    for log_file in log_files:
        try:
            # Need to add an entry that corresponds to the original log file and the original index.
            log_df = pd.read_csv(log_file)
            log_df["log"] = [log_file]*len(log_df)
            log_df["id"] = list(log_df.index)
        except pd.errors.EmptyDataError:
            continue
        log_dfs.append(log_df[equil::offset])

    # Return np.nan if there are no log files or the length of the merged log file is 0.
    if len(log_dfs) == 0:
        field_stat_df = pd.DataFrame(columns=["field", "stat", "value", "log", "id"])
        field_stat_df.loc[0] = [field, stat, np.nan, np.nan, np.nan]
    else:
        merge_log_df = pd.concat(log_dfs)
        field_stat_df = get_stat_from_log_df(
                log_df=merge_log_df,
                field=field,
                stat=stat,
                N=N
            )

    return log_files, field_stat_df


"""
get_stat_df function takes in a set of output directories (presumably the set of output directories for a job) and returns a dataframe containing the stat for each field for each log file for all stat, field pair. Having this dataframe is important for sample volume benchmarks.

params
log_file_groups: the set of groups where each group contains one or more log files (the rows of the stat_df). Each stat, field pair is computed for each group of log files.
fields: the fields (eg, rmsd) from each log file to extract a statistic for.
stats: the statistic (eg, mean) to compute for each field (field[i]_statistic[i] is a column).
N: the number of entries to return for each field, stat pair. This only applies to min and max.
offset: only look return stats for log entries with offset frequency. This is useful when we only want entries that have a corresponding pdb file.
equil: the number of frames to discard before computing a statistic.
test: additional param to get the stat_df without multiprocessing.

returns
stat_df: pandas dataframe, dataframe with each log file group as entires and the requested stats as columns. The type and number of returned columns depends on the stat type as well as other parameters.
"""
def get_stat_df(
        log_file_groups,
        fields,
        stats,
        N,
        offset,
        equil,
        test=False
):
    if len(fields) != len(stats):
        raise RuntimeError("The number of fields ({}) does not match the number of stats ({})".format(len(fields), len(stats)))

    n_fields = len(fields)

    # Construct the columns for the stat_df.
    columns = list()
    for i in range(n_fields):
        if stats[i] in ["min", "max"]:
            for j in range(N):
                column = "{}_{}_{}".format(fields[i], stats[i], j)
                columns.append(column)
                columns.append("{}_log".format(column))
                columns.append("{}_id".format(column))
        else:
            column = "{}_{}".format(fields[i], stats[i])
            columns.append(column)

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
            pool_param["N"] = N
            pool_param["offset"] = offset
            pool_params.append(pool_param)

    # print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    if test:
        pool_results = list()
        for pool_param in pool_params:
            pool_results.append(pool_get_stat_from_log_files(pool_param))
    else:
        pool_results = pool_obj.imap(
            pool_get_stat_from_log_files,
            pool_params
        )

    # Iterate through the pool results and fill in the stat_df.
    # print(log_stat_df.head())
    for log_files, field_stat_df in pool_results:
        # print(field_stat_df.head())
        field = field_stat_df.iloc[0]["field"]
        stat = field_stat_df.iloc[0]["stat"]

        entry = str(log_files)

        for i in range(len(field_stat_df)):
            if stat in ["min", "max"]:
                col = "{}_{}_{}".format(field, stat, i)
            else:
                col = "{}_{}".format(field, stat, i)

            val = field_stat_df.iloc[i]["value"]

            log_stat_df.loc[entry, col] = val

            log_file_col = "{}_log".format(col)
            id_col = "{}_id".format(col)
            if log_file_col in log_stat_df.columns:
                log_file = field_stat_df.iloc[i]["log"]
                entry_id = field_stat_df.iloc[i]["id"]

                log_stat_df.loc[entry, log_file_col] = log_file
                log_stat_df.loc[entry, id_col] = entry_id

    pool_obj.close()
    return log_stat_df
