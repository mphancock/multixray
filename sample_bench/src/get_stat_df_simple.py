from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing


def pool_get_stat_info_df(
        params_dict
):
    log_files = params_dict["log_files"]
    equil = params_dict["equil"]
    field = params_dict["field"]
    bonus_fields = params_dict["bonus_fields"]

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

    N = params_dict["N"]
    pdb_only = params_dict["pdb_only"]
    max_rmsd = params_dict["max_rmsd"]

    # Merge the log files into a single dataframe.
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
        check_fields = list()
        check_fields.extend(bonus_fields)
        check_fields.extend(fields)

        for i in range(len(check_fields)):
            tmp_field = check_fields[i]
            if tmp_field not in log_df.columns:
                print("{} missing field: {}".format(log_file, tmp_field))
                invalid = True

        if invalid:
            continue

        if pdb_only:
            log_df = log_df[~log_df['pdb'].isna()]

        log_df = log_df.iloc[equil:]

        if max_rmsd:
            log_df = log_df[log_df["rmsd_0"] <= max_rmsd]

        log_dfs.append(log_df)

    columns = [field]
    columns.extend(params_dict["bonus_fields"])
    columns = list(set(columns))
    merge_log_df = pd.DataFrame(columns=columns)

    if len(log_dfs) > 0:
        merge_log_df = pd.concat(log_dfs)

        # Merge back + fields and bonus fields.
        if len(fields) > 1:
            merge_log_df[field] = merge_log_df[fields].sum(axis=1)

        for bonus_field in params_dict["bonus_fields"]:
            if "+" in bonus_field:
                merge_log_df[bonus_field] = merge_log_df[bonus_field.split("+")].sum(axis=1)

        if len(merge_log_df) > (equil+N):
            stat_df = merge_log_df.nsmallest(N, field)
            stat_df = stat_df[columns]
        else:
            print("Not enough entries to compute stat: {}".format(log_files))
            stat_df = merge_log_df
    else:
        print("No valid log dfs")
        stat_df = merge_log_df

    return stat_df


def get_stat_df(
        log_files,
        field,
        N,
        bonus_fields,
        equil,
        pdb_only,
        max_rmsd=None
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
        pool_params.append(pool_param)

    # Collect all of the field stat dfs and log file groups that are returned by get_stat_from_log_files. The field stat dfs will be used to populate the final stat df. The format of the field stat dfs has the following columns: field, stat, value, log, id. The df will only have more than one row if the stat is either "max" or "min" and N>1.
    pool_results = list()
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
    job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/166_N1/6")
    # log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("*/")]
    log_files = [Path(job_dir, "output_0/log.csv")]
    print(log_files)

    stat_df = get_stat_df(
        log_files=log_files,
        field="xray_0+xray_1",
        N=1,
        bonus_fields=["pdb", "ff"],
        equil=1,
        pdb_only=True
    )

    print(stat_df.columns)
    print(stat_df.head())

    # print(stat_df.columns)
    # print(stat_df.iloc[0, 0])
    # print(stat_df.iloc[0, 1])
    # print(stat_df.iloc[0, 2])