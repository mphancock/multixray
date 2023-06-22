from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing

"""
[Path(/abc/def), Path(/ghi/jkl)] -> "/abc/def,/ghi/jkl"
"""
def get_log_file_group_str(
        log_file_group
):
    log_file_group_str = "{}".format(log_file_group[0])
    for i in range(1, len(log_file_group)):
        log_file = log_file_group[i]
        log_file_group_str = log_file_group_str + ",{}".format(str(log_file))

    return log_file_group_str

"""
pool_get_stat_info_df returns the stat_info_df from a single log_df. The stat_inf_df is an inversion of the stat_df (kind of) such that the stat_df has the log_file_groups as rows while the stat_info_df has the statistics as rows. If the statistic is in mean, std, or var; then there is only 1 row and the only field is the statistic value. If the statistic is min or max, then there are N rows, and the fields are the field + bonus field values from the min/maximizing log_df entries.

**********
Parameters
    params_dict: parameter dictionary containing parameters for reading and analyzing the log files.

**********
Returns
    log_files: list of Path objects, the log files that were read to compute the statistics.

    field_stat_df: pandas dataframe, dataframe containing the stat value(s) plus entry information. The dataframe will only be greater than a single row for min and max with N>1. The columns are [main_field, *bonus_fields] and the rows are [0, 1, ..., N-1].



"""
def pool_get_stat_info_df(
        params_dict
):
    log_files = params_dict["log_files"]
    equil = params_dict["equil"]
    field = params_dict["field"]
    bonus_fields = params_dict["bonus_fields"]
    stat = params_dict["stat"]
    N = params_dict["N"]
    offset = params_dict["offset"]
    max_frame = params_dict["max_frame"]

    # Merge the log files into a single dataframe.
    log_dfs = list()
    for log_file in log_files:
        try:
            log_df = pd.read_csv(log_file)
        except pd.errors.EmptyDataError:
            print("Skipped empty log_file: {}".format(log_file))
            continue

        pdb_files = list()
        for i in range(len(log_df)):
            step = log_df["step"].iloc[i]
            if step % 10 == 0:
                pdb_files.append(str(Path(log_file.parents[0], "pdbs", "{}.pdb".format(step // 10))))
            else:
                pdb_files.append(np.nan)
        log_df["pdb"] = pdb_files

        if max_frame is not None:
            log_df_subset = log_df[equil:max_frame:offset]
        else:
            log_df_subset = log_df[equil::offset]
        log_dfs.append(log_df_subset)


    columns = [field]
    columns.extend(bonus_fields)
    stat_info_df = pd.DataFrame(columns=columns)

    if len(log_dfs) > 0:
        merge_log_df = pd.concat(log_dfs)
        if len(merge_log_df) > 0:
            # Check that all the fields, bonus fields, and stats are valid.
            if stat not in ["min", "max", "mean", "std", "var"]:
                return RuntimeError("Stat {} not valid".format(stat))
            for column in columns:
                if column not in merge_log_df.columns:
                    return RuntimeError("Column {} not in merge_log_df".format(column))

            # Compute the stat. Return np.nan if there are no log files or the length of the merged log file is 0.
            if stat in ["min", "max"]:
                if stat == "min":
                    n_stat_df = merge_log_df.nsmallest(N, field)
                else:
                    n_stat_df = merge_log_df.nlargest(N, field)

                stat_info_df = n_stat_df[columns]
            else:
                if stat == "mean":
                    stat_val = merge_log_df[field].mean()
                elif stat == "std":
                    stat_val = merge_log_df[field].std()
                elif stat == "var":
                    stat_val = merge_log_df[field].var()

                stat_info_df.loc[0] = [stat_val]
        else:
            stat_info_df.loc[0] = [np.nan]*len(stat_info_df.columns)
    else:
        stat_info_df.loc[0] = [np.nan]*len(stat_info_df.columns)

    return log_files, stat_info_df


"""
get_stat_df function takes in a set of output directories (presumably the set of output directories for a job) and returns a dataframe containing the stat for each field for each log file for all stat, field pair. Having this dataframe is important for sample volume benchmarks.

*********
params
    log_file_groups: the set of groups where each group contains one or more log files (the rows of the stat_df). Each stat, field pair is computed for each group of log files.

    main_field: the field to compute the statistic for.

    stats: the statistic (eg, mean) to compute.

    N: the number of entries to return for each field, stat pair. This only applies to min and max.

    offset: only look return stats for log entries with offset frequency. This is useful when we only want entries that have a corresponding pdb file.

    equil: the number of frames to discard before computing a statistic.

    max_frame: the number of frames to consider when computing a statistic. This is useful when there are more frames than there are corresponding pdb files and we do not want to consider surplus frames.

    test: additional param to get the stat_df without multiprocessing.

*********
Return
    stat_df: pandas dataframe containing the requested statistics. The rows are the log_file_groups. If the main_field is an aggregrate, then the columns are [field_value]. If the main_field is min/max then the columns are [field_value_i, field_value_i_*bonus_columns] for 0<i<N.

"""
def get_stat_df(
        log_file_groups,
        main_field,
        main_stat,
        N=1,
        bonus_fields=[],
        offset=1,
        equil=0,
        max_frames=None,
        test=False
):
    # Construct the stat_df.
    columns = list()

    if main_stat in ["min", "max"]:
        for i in range(N):
            col = "{}_{}_{}".format(main_field, main_stat, i)
            columns.append(col)

            for bonus_field in bonus_fields:
                columns.append("{}_{}".format(col, bonus_field))
    else:
        columns.append("{}_{}".format(main_field, main_stat))

    indices = [get_log_file_group_str(group) for group in log_file_groups]
    log_stat_df = pd.DataFrame(index=indices, columns=columns)

    # If there is only a single file group and a single min or max stat, then partition the single file group into batches to distribute computation.
    fast = False
    if len(log_file_groups) == 1 and main_stat in ["min", "max"]:
        fast = True

        if multiprocessing.cpu_count() > len(log_file_groups[0]):
            n_batches = len(log_file_groups[0])
        else:
            n_batches = multiprocessing.cpu_count()

        log_file_groups_tmp = np.array_split(log_file_groups[0], n_batches)
    else:
        log_file_groups_tmp = log_file_groups

    # Iterate through each group of log files to distribute the calculation of a each field, stat pair.
    pool_params = list()
    # for log_file_group in log_file_groups_tmp:
    for i in range(len(log_file_groups_tmp)):
        log_file_group = log_file_groups_tmp[i]

        # This doesn't work with fast.
        if max_frames and not fast:
            max_frame = max_frames[i]
        else:
            max_frame = None

        # For each field, stat pair, create a pool_param.
        pool_param = dict()
        pool_param["log_files"] = log_file_group
        pool_param["equil"] = equil
        pool_param["field"] = main_field
        pool_param["stat"] = main_stat
        pool_param["bonus_fields"] = bonus_fields
        pool_param["N"] = N
        pool_param["offset"] = offset
        pool_param["max_frame"] = max_frame
        pool_params.append(pool_param)

    # Collect all of the field stat dfs and log file groups that are returned by get_stat_from_log_files. The field stat dfs will be used to populate the final stat df. The format of the field stat dfs has the following columns: field, stat, value, log, id. The df will only have more than one row if the stat is either "max" or "min" and N>1.
    results = list()
    if test:
        for pool_param in pool_params:
            results.append(pool_get_stat_info_df(pool_param))
    else:
        pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
        pool_results = pool_obj.imap(pool_get_stat_info_df, pool_params)

        for pool_result in pool_results:
            if isinstance(pool_result, Exception):
                pool_obj.close()
                raise pool_result
            else:
                results.append(pool_result)

    # If the fast flag is turned on, then the field_stat_dfs, and log_files need to be merged.
    if fast:
        batch_stat_info_dfs = list()
        for result in results:
            batch_log_files, batch_stat_info_df = result
            batch_stat_info_dfs.append(batch_stat_info_df)

        merge_stat_info_df = pd.concat(batch_stat_info_dfs)
        if main_stat == "min":
            stat_info_df = merge_stat_info_df.nsmallest(N, main_field)
        elif main_stat == "max":
            stat_info_df = merge_stat_info_df.nlargest(N, main_field)
        else:
            raise RuntimeError("Something has gone wrong")

        results = [(log_file_groups[0], stat_info_df)]

    for log_file_group, stat_info_df in results:
        entry = get_log_file_group_str(log_file_group)

        if main_stat in ["mean", "std", "var"]:
            col = "{}_{}".format(main_field, main_stat)
            log_stat_df.loc[entry,col] = stat_info_df[main_field].iloc[0]
        else:
            for i in range(N):
                col = "{}_{}_{}".format(main_field, main_stat, i)
                log_stat_df.loc[entry,col] = stat_info_df[main_field].iloc[i]

                for j in range(len(bonus_fields)):
                    bonus_col = col+"_{}".format(bonus_fields[j])
                    log_stat_df.loc[entry,bonus_col] = stat_info_df[bonus_fields[j]].iloc[i]

    if not test:
        pool_obj.close()

    return log_stat_df