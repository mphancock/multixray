from pathlib import Path
import random
import pandas as pd
random.seed(0)
import IMP
import IMP.atom
import sys
import multiprocessing
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/src")))
import weights


def get_random_sample_df(
    log_dfs,
    N,
    equil,
    range,
    field
):
    pdb_log_dfs = list()

    for log_df in log_dfs:
        pdb_log_df = log_df[log_df['pdb'].notna()]

        # remove the first structure.
        pdb_log_df = pdb_log_df.iloc[equil:]
        pdb_log_df = pdb_log_df[(pdb_log_df[field] >= range[0]) & (pdb_log_df[field] <= range[1])]
        pdb_log_dfs.append(pdb_log_df)

    merge_log_df = pd.concat(pdb_log_dfs)
    print(len(merge_log_df))

    if len(merge_log_df) < N:
        sample_log_df = merge_log_df
    else:
        sample_log_df = merge_log_df.sample(n=N)

    return sample_log_df



def get_valid_output_dirs(
        job_dirs
):
    output_dirs = list()
    for job_dir in job_dirs:
        for out_dir in job_dir.glob("output_*"):
            print(out_dir)
            if Path(out_dir, "pdbs").exists() and Path(out_dir, "log.csv").exists():
                output_dirs.append(out_dir)

    return output_dirs


def read_log_file(
    log_file
):
    # Log file may be empty.
    try:
        log_df = pd.read_csv(log_file, index_col=0)
    except pd.errors.EmptyDataError as e:
        return e
    except FileNotFoundError as e:
        return e

    return log_file, log_df


def get_all_log_dfs(
        log_files
):
    log_dfs = list()
    pool_obj = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool_results = pool_obj.imap(
        read_log_file,
        log_files
    )

    for result in pool_results:
        if isinstance(result, Exception):
            continue
        else:
            log_file, log_df = result
            log_dfs.append(log_df)

    return log_dfs


if __name__ == "__main__":
    exp_name = "268_decoys_1_state"

    rmsd_ranges = [[0,0.25],[.25,1.0]]
    n_decoys = [750, 250]

    job_dir = Path(Path.home(), "xray/score_bench/data", exp_name)
    job_dir.mkdir(exist_ok=True)

    n_state = 1
    n_cond = 1
    decoy_meta_file = Path(job_dir, "rand1000.csv")

    sample_job_dirs = [Path("/wynton/group/sali/mhancock/xray/sample_bench/out/{}/0".format(exp_name)),Path("/wynton/group/sali/mhancock/xray/sample_bench/out/{}/1".format(exp_name))]

    out_dirs = get_valid_output_dirs(job_dirs=sample_job_dirs)
    log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]

    # We should read the log_df here so that it isn't done for every RMSD range.
    log_dfs = get_all_log_dfs(log_files=log_files)

    sample_dfs = list()
    for i in range(len(rmsd_ranges)):
        print(rmsd_ranges[i])
        # Get and save the dataframe containing all the decoy entries (each decoy entry containing 1 or more random pdb files and an equal number of weights).
        sample_df = get_random_sample_df(
            log_dfs=log_dfs,
            N=n_decoys[i],
            equil=1,
            range=rmsd_ranges[i],
            field="rmsd"
        )
        sample_dfs.append(sample_df)

    decoy_df = pd.concat(sample_dfs)

    cols = ["pdb", "ff", "rmsd"]
    for state in range(n_state):
        for cond in range(n_cond):
            cols.append("w_{}_{}".format(state, cond))
    decoy_df = decoy_df[cols]

    # decoy_df.drop(columns=["time", "step", "copy"], inplace=True)
    decoy_df.reset_index(drop=True, inplace=True)
    decoy_df.to_csv(decoy_meta_file)
