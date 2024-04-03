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

    if "+" in field:
        fields = field.split("+")
    else:
        fields = [field]

    # Check the bonus field for +.
    bonus_fields = list()
    for bonus_field in params_dict["bonus_fields"]:
        if "+" in bonus_field:
            bonus_fields.extend(bonus_field.split("+"))
        else:
            bonus_fields.append(bonus_field)

    stat = params_dict["stat"]
    N = params_dict["N"]
    pdb_only = params_dict["pdb_only"]
    rmsd_filter = params_dict["rmsd_filter"]

    # Merge the log files into a single dataframe.
    log_dfs = list()
    for log_file in log_files:
        try:
            log_df = pd.read_csv(log_file)
        except pd.errors.EmptyDataError:
            print("Skipped empty log_file: {}".format(log_file))
            continue
        except FileNotFoundError:
            print("Skipped missing log_file: {}".format(log_file))
            continue

        # Create a set of columns to check for existence
        check_cols = ["step", "time"]
        check_cols.extend(bonus_fields)
        check_cols.extend(fields)

        # Check if all columns in check_cols exist in the DataFrame. This guaruntees that merge_log_df will contain all of the necessary columns.
        df_contains_fields = True
        for col in check_cols:
            if col not in log_df.columns:
                df_contains_fields = False
                print("Skipped {} missing fields: {}".format(log_file, col))
                break

        if not df_contains_fields:
            continue

        if pdb_only:
            log_sel_df = log_df[~log_df['pdb'].isna()].iloc[equil:]
        else:
            log_sel_df = log_df.iloc[equil:]

        if rmsd_filter is not None:
            log_sel_df = log_sel_df[log_sel_df["rmsd_0"] <= rmsd_filter]

        log_dfs.append(log_sel_df)

    columns = [field]
    columns.extend(params_dict["bonus_fields"])
    stat_info_df = pd.DataFrame(columns=columns)

    if len(log_dfs) > 0:
        merge_log_df = pd.concat(log_dfs)

        # Merge back + fields and bonus fields.
        if len(fields) > 1:
            merge_log_df[field] = merge_log_df[fields].sum(axis=1)

        for bonus_field in params_dict["bonus_fields"]:
            if "+" in bonus_field:
                merge_log_df[bonus_field] = merge_log_df[bonus_field.split("+")].sum(axis=1)

        if len(merge_log_df) > (equil+N):
            # Check that all the fields, bonus fields, and stats are valid.
            if stat not in ["min", "max", "mean", "std", "var"]:
                return RuntimeError("Stat {} not valid".format(stat))

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
            print("Not enough entries to compute stat: {}".format(log_files))

    return log_files, stat_info_df


"""
get_stat_df function takes in a set of output directories (presumably the set of output directories for a job) and returns a dataframe containing the stat for each field for each log file for all stat, field pair. Having this dataframe is important for sample volume benchmarks.

With the updates to the logging we now know both which log entries have associated pdb files, as well as the # of pdb files being <= what is recorded in the log.

*********
params
    log_file_groups: the set of groups where each group contains one or more log files (the rows of the stat_df). Each stat, field pair is computed for each group of log files.

    main_field: the field to compute the statistic for.

    stats: the statistic (eg, mean) to compute.

    N: the number of entries to return for each field, stat pair. This only applies to min and max.

    equil: the number of frames to discard before computing a statistic.

    pdb_only: only compute the statistic for log entries that have a pdb file.

    rmsd_filter: only compute the statistic for log entries that have an rmsd less than or equal to the rmsd_filter. This is useful for computing statistics from log files with unphysiological structures.

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
        equil=0,
        pdb_only=False,
        rmsd_filter=None,
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
    for i in range(len(log_file_groups_tmp)):
        log_file_group = log_file_groups_tmp[i]

        # For each field, stat pair, create a pool_param.
        pool_param = dict()
        pool_param["log_files"] = log_file_group
        pool_param["equil"] = equil
        pool_param["field"] = main_field
        pool_param["stat"] = main_stat
        pool_param["bonus_fields"] = bonus_fields
        pool_param["N"] = N
        pool_param["pdb_only"] = pdb_only
        pool_param["rmsd_filter"] = rmsd_filter
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
                log_file_group, stat_info_df = pool_result

                # There are instances where the length of stat_info_df will be 0 (eg, if the log files are empty)
                if len(stat_info_df) > 0:
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


if __name__ == "__main__":
    stat_df = get_stat_df(
        log_file_groups=[[Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/158_N8_J3/640510/output_9/log.csv")]],
        main_field="xray_0+xray_1+xray_2",
        main_stat="min",
        bonus_fields=["pdb", "ff"],
        pdb_only=True,
        equil=1
    )

    print(stat_df.head())

    # print(stat_df.columns)
    # print(stat_df.iloc[0, 0])
    # print(stat_df.iloc[0, 1])
    # print(stat_df.iloc[0, 2])