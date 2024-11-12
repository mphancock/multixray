from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing


def pool_get_stat_info_df(
        params_dict
):
    log_files = params_dict["log_files"]
    equil = params_dict["equil"]

    ## field and bonus field are the orignal field and bonus field passed in the parameters dict and may include sum fields (eg, xray_0+xray_1)
    field = params_dict["field"]
    bonus_fields = params_dict["bonus_fields"]

    if "+" in field:
        tmp_fields = field.split("+")
    else:
        tmp_fields = [field]

    ## tmp fields and tmp bonus fields are the fields that are actually present in the log files (eg, xray_0, xray_1)
    tmp_bonus_fields = list()
    for bonus_field in bonus_fields:
        if "+" in bonus_field:
            tmp_bonus_fields.extend(bonus_field.split("+"))
        else:
            tmp_bonus_fields.append(bonus_field)

    N = params_dict["N"]
    pdb_only = params_dict["pdb_only"]
    max_rmsd = params_dict["max_rmsd"]
    max_ff = params_dict["max_ff"]

    # merge the log files into a single dataframe.
    log_dfs = []
    for log_file in log_files:
        try:
            log_df = pd.read_csv(log_file)
        except pd.errors.EmptyDataError:
            print("Skipped empty log_file: {}".format(log_file))
            continue
        except FileNotFoundError:
            print("Skipped missing log_file: {}".format(log_file))
            continue

        invalid = False
        check_fields = tmp_fields.copy()
        check_fields.extend(tmp_bonus_fields)

        for i in range(len(check_fields)):
            tmp_field = check_fields[i]
            if tmp_field not in log_df.columns:
                print("{} missing field: {}".format(log_file, tmp_field))
                invalid = True

        if invalid:
            continue

        log_df = log_df.iloc[equil:]

        if pdb_only:
            log_df = log_df[~log_df['pdb'].isna()]

        if max_rmsd:
            log_df = log_df[log_df["rmsd_0"] <= max_rmsd]

        if max_ff:
            log_df = log_df[log_df["ff"] <= max_ff]

        log_dfs.append(log_df)

    merge_log_df = pd.DataFrame()
    if len(log_dfs) > 0:
        merge_log_df = pd.concat(log_dfs)

        # compute the sum fields
        if len(tmp_fields) > 1:
            merge_log_df[field] = merge_log_df[tmp_fields].sum(axis=1)

        for bonus_field in bonus_fields:
            if "+" in bonus_field:
                merge_log_df[bonus_field] = merge_log_df[bonus_field.split("+")].sum(axis=1)

        if len(merge_log_df) >= (equil+N):
            stat_df = merge_log_df.nsmallest(N, field)
        else:
            print("Not enough entries to compute stat: {}".format(log_files))
            stat_df = merge_log_df
    else:
        print("No valid log dfs")
        stat_df = merge_log_df

    ## columns to ask for in the final stat_df
    columns = [field]
    columns.extend(bonus_fields)
    columns = list(set(columns))

    stat_df = stat_df[columns]

    return stat_df


def get_stat_df(
        log_files,
        field,
        N,
        bonus_fields,
        equil,
        pdb_only,
        max_rmsd=None,
        max_ff=None
):
    # Construct the stat_df.
    columns = [field]
    columns.extend(bonus_fields)

    fast = False

    n_cpu = multiprocessing.cpu_count()
    if len(log_files) > n_cpu:
        n_batch = n_cpu
    else:
        n_batch = len(log_files)

    log_file_groups = np.array_split(log_files, n_batch)

    # Iterate through each group of log files to distribute the calculations.
    pool_params = list()
    for log_file_group in log_file_groups:
        pool_param = dict()
        pool_param["log_files"] = log_file_group
        pool_param["equil"] = equil
        pool_param["field"] = field
        pool_param["N"] = N
        pool_param["bonus_fields"] = bonus_fields
        pool_param["pdb_only"] = pdb_only
        pool_param["max_rmsd"] = max_rmsd
        pool_param["max_ff"] = max_ff
        pool_params.append(pool_param)

    # Collect all of the field stat dfs and log file groups that are returned by get_stat_from_log_files. The field stat dfs will be used to populate the final stat df. The format of the field stat dfs has the following columns: field, stat, value, log, id. The df will only have more than one row if the stat is either "max" or "min" and N>1.
    # pool_results = list()
    # for pool_param in pool_params:
    #     pool_results.append(pool_get_stat_info_df(pool_param))

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(pool_get_stat_info_df, pool_params)

    subset_stat_dfs = list()
    for pool_result in pool_results:
        subset_stat_dfs.append(pool_result)

    pool_obj.close()

    merge_stat_df = pd.concat(subset_stat_dfs)
    stat_df = merge_stat_df.nsmallest(N, field)

    stat_df = stat_df.reset_index()

    return stat_df


if __name__ == "__main__":
    job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/267_full_ref/0")
    log_files = [Path(job_dir, "output_0/log.csv")]
    print(log_files)

    stat_df = get_stat_df(
        log_files=log_files,
        field="xray_native_0+xray_native_1",
        N=1,
        bonus_fields=["pdb", "ff", "rmsd_native_0+rmsd_native_1"],
        equil=0,
        pdb_only=True
    )

    print(stat_df.columns)
    print(stat_df.head())