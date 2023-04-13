from pathlib import Path
import argparse
import sys
import pandas as pd
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import log_file_analysis


""" get_stat_df_from_logs function takes in a set of output directories (presumably the set of output directories for a job) and returns a dataframe containing the minimum value for each specified field for each log file. Having this dataframe is important for sample volume benchmarks.

log_files: the set of log files to compute statistics from (rows).  
fields: the fields (eg, rmsd) from each log file to extract a statistic for. 
statistics: the statistic (eg, mean) to compute for each field (field[i]_statistic[i] is a column).  
equil: the number of frames to discard before computing a statistic. 
include_file: whether to include an additional column for each field, statistic pair which records the relevant file. This can only be used for min and max statistics. 
"""
def get_stat_df_from_logs(
        log_files,
        fields,
        stats,
        equil,
        include_id
):
    pool_params = list()
    n_fields = len(fields)
    valid_log_files = list()
    for log_file in log_files:
        try:
            log_df = pd.read_csv(log_file)
        except pd.errors.EmptyDataError:
            continue

        if len(log_df) < equil:
            continue

        valid_log_files.append(log_file)
        log_df_equil = log_df.iloc[equil:]

        for i in range(n_fields):
            pool_param = dict()
            pool_param["log_file"] = log_file

            pool_param["log_df"] = log_df_equil
            pool_param["field"] = fields[i]
            pool_param["stat"] = stats[i]
            pool_params.append(pool_param)
    n_log_files = len(valid_log_files)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        log_file_analysis.pool_get_stat_from_log_df,
        pool_params
    )

    columns = ["log_file"]
    for i in range(n_fields):
        columns.append("{}_{}".format(fields[i], stats[i]))
        if stats[i] in ["min", "max"] and include_id:
            columns.append("{}_{}_id".format(fields[i], stats[i]))

    log_stat_df = pd.DataFrame(index=list(range(n_log_files)), columns=columns)

    for i in range(n_log_files):
        log_stat_df.loc[i, "log_file"] = valid_log_files[i]

    for pool_result in pool_results:
        if not pool_result:
            continue
        else:
            log_file, field, stat, log_stat, stat_id = pool_result

        log_file_id = list(log_stat_df[log_stat_df['log_file'] == log_file].index)[0]
        log_stat_df.iloc[log_file_id, log_stat_df.columns.get_loc("{}_{}".format(field, stat))] = log_stat
        if stat_id and include_id:
            log_stat_df.iloc[log_file_id, log_stat_df.columns.get_loc("{}_{}_id".format(field, stat))] = stat_id

    return log_stat_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--job_name")
    parser.add_argument("--job_id", type=int)
    args = parser.parse_args()
    print(args.job_name)
    print(args.job_id)

    job_name = args.job_name
    exp_id = job_name.split("_")[0]
    job_id = args.job_id

    xray_dir = Path("/wynton/group/sali/mhancock/xray")
    out_dirs = list(Path(xray_dir, "sample_bench/out/3ca7/single_md/{}/{}".format(job_name, job_id)).glob("output*"))
    fields = ["tot", "rmsd"]

    # First, construct the lookup table that contains for each log_file entry, the minimum value for that MD log for each of the requested fields.
    log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]
    stat_df = get_stat_df_from_logs(
        log_files=log_files,
        fields=["tot", "tot", "rmsd", "rmsd"],
        stats=["min", "mean", "min", "mean"],
        equil=1000,
        include_id=True
    )
    print(stat_df.head())

    # Create an additional column that is the rmsd of the minimum scoring structure.
    min_score_rmsds = list()
    min_score_rmsd_aligns = list()
    for i in range(len(stat_df)):
        log_file = stat_df.loc[i, "log_file"]
        min_score_id = stat_df.loc[i, "tot_min_id"]

        log_df = pd.read_csv(log_file, index_col=0)
        min_score_rmsd = log_df.loc[min_score_id, "rmsd"]
        min_score_rmsd_align = log_df.loc[min_score_id, "rmsd_align"]
        min_score_rmsds.append(min_score_rmsd)
        min_score_rmsd_aligns.append(min_score_rmsd_align)

    stat_df["tot_min_rmsd"] = min_score_rmsds
    stat_df["tot_min_rmsd_align"] = min_score_rmsd_aligns

    stat_df.to_csv(Path(Path.home(), "xray/sample_bench/data/stat_df/{}.csv".format(exp_id)))
