from pathlib import Path
import random
import pandas as pd
random.seed(0)
import IMP
import IMP.atom
import sys
import multiprocessing
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
import weights


def get_random_sample_df(
    log_dfs,
    N,
    equil,
    rmsd_range
):
    pdb_log_dfs = list()

    for log_df in log_dfs:
        pdb_log_df = log_df[log_df['pdb'].notna()]

        # remove the first structure.
        pdb_log_df = pdb_log_df.iloc[1:]

        if equil:
            pdb_log_df = pdb_log_df.iloc[int(equil*len(pdb_log_df)):]

        if rmsd_range:
            pdb_log_df = pdb_log_df[(pdb_log_df['rmsd_avg_0'] >= rmsd_range[0]) & (pdb_log_df['rmsd_avg_0'] <= rmsd_range[1])]

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
    target = "3ca7"
    job_name = "152_native_1_cif"

    n_decoys = 1000
    rmsd_ranges = list()

    rmsd_ranges.append([0,1])
    # rmsd_ranges.append([.5,1])
    job_dir = Path(Path.home(), "xray/score_bench/data", target, job_name)
    job_dir.mkdir(exist_ok=True)

    n_cif = 1
    n_state = 2
    fields = ["pdb", "ff"]

    # ###
    # fields.append("rmsd_avg_0")

    # decoy_meta_file = Path(job_dir, "rand1000.csv")

    # sample_job_dirs = list()
    # sample_job_dirs.append(Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/54_1000/4474326"))
    # out_dirs = get_valid_output_dirs(job_dirs=sample_job_dirs)
    # log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]

    # # We should read the log_df here so that it isn't done for every RMSD range.
    # log_dfs = get_all_log_dfs(log_files=log_files)

    # sample_dfs = list()
    # for i in range(len(rmsd_ranges)):
    #     print(rmsd_ranges[i])
    #     # Get and save the dataframe containing all the decoy entries (each decoy entry containing 1 or more random pdb files and an equal number of weights).
    #     sample_df = get_random_sample_df(
    #         log_dfs=log_dfs,
    #         N=n_decoys,
    #         equil=0,
    #         rmsd_range=rmsd_ranges[i]
    #     )

    #     sample_df = sample_df[fields]
    #     sample_dfs.append(sample_df)

    # decoy_df = pd.concat(sample_dfs)
    # decoy_df.reset_index(drop=True, inplace=True)
    # decoy_df.to_csv(decoy_meta_file)

    ###

    for i in range(n_cif):
        fields.append("xray_{}".format(i))
        fields.append("r_free_{}".format(i))
        fields.append("rmsd_avg_{}".format(i))

        if n_state > 1:
            for j in range(n_state):
                fields.append("weight_{}_{}".format(i, j))

    for job_id in range(10):
        decoy_id = 0
        decoy_name = str(job_id)

        decoy_pdb_dir = Path("/wynton/group/sali/mhancock/xray/decoys/data", target, job_name, decoy_name)
        decoy_pdb_dir.mkdir(exist_ok=True, parents=True)

        decoy_meta_file = Path(job_dir, "{}.csv".format(decoy_name))

        sample_job_dirs = list()
        sample_job_dirs.append(Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/{}/{}".format(job_name, job_id)))
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
                N=n_decoys,
                equil=.5,
                rmsd_range=rmsd_ranges[i]
            )

            sample_df = sample_df[fields]
            sample_dfs.append(sample_df)

        decoy_df = pd.concat(sample_dfs)
        decoy_df.reset_index(drop=True, inplace=True)
        decoy_df.to_csv(decoy_meta_file)