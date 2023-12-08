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
import merge_pdbs


def get_random_sample_df(
    log_dfs,
    N,
    equil,
    rmsd_range
):
    pdb_log_dfs = list()

    for log_df in log_dfs:
        pdb_log_df = log_df[log_df['pdb'].notna()].iloc[equil:]
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


def get_occs(
        pdb_file
):
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    occs = list()
    for h in hs:
        pids = IMP.atom.Selection(h).get_selected_particles()
        occ = IMP.atom.Atom(m, pids[0]).get_occupancy()
        occs.append(occ)

    return occs


def get_n_states(
        pdb_file
):
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    return len(hs)


def sample_occs(
        n_decoys,
        n_state,
        occ_means,
        rmsd_range
):
    all_occs = list()
    for i in range(n_decoys):
        occs = list()
        for j in range(n_state):
            rand_num = np.random.normal(occ_means[j], 1/10*rmsd_range[0], n_state)[0]

            if rand_num < 0:
                rand_num = 0

            occs.append(rand_num)

        all_occs.append(occs)

    # print(all_occs)

    # Normalize.
    all_occs_norm = list()
    for occs in all_occs:
        occs_norm = [occ/sum(occs) for occ in occs]
        all_occs_norm.append(occs_norm)

    return all_occs_norm


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


def get_occs_from_pdb_files(
        pdb_files
):
    pool_obj = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool_results = pool_obj.imap(
        get_occs,
        pdb_files
    )

    occs = list()
    for result in pool_results:
        occs.append(result)

    return occs



if __name__ == "__main__":
    target = "3ca7"
    job_name = "100_natives_4x"

    n_struct = 1
    n_decoys = 1000
    rmsd_ranges = list()

    rmsd_ranges.append([0, 3])
    # for i in range(10):
    #     rmsd_ranges.append([i*.1, (i+1)*.1])

    decoy_meta_dir = Path(Path.home(), "xray/decoys/data", target, job_name)
    decoy_meta_dir.mkdir(exist_ok=True, parents=True)

    n_cif = 1
    n_state = 4
    fields = ["pdb"]

    for i in range(n_cif):
        fields.append("xray_{}".format(i))
        fields.append("r_free_{}".format(i))
        fields.append("rmsd_avg_{}".format(i))
        for j in range(n_state):
            fields.append("weight_{}_{}".format(i, j))

    for job_id in range(10):
        decoy_id = 0
        decoy_name = str(job_id)

        decoy_pdb_dir = Path("/wynton/group/sali/mhancock/xray/decoys/data", target, job_name, decoy_name)
        decoy_pdb_dir.mkdir(exist_ok=True, parents=True)

        decoy_meta_file = Path(decoy_meta_dir, "{}.csv".format(decoy_name))

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
                equil=1,
                rmsd_range=rmsd_ranges[i]
            )

            sample_df = sample_df[fields]
            sample_dfs.append(sample_df)

        decoy_df = pd.concat(sample_dfs)
        decoy_df.reset_index(drop=True, inplace=True)
        decoy_df.to_csv(decoy_meta_file)