from pathlib import Path
import random
import pandas as pd
random.seed(0)
import IMP
import IMP.atom
import sys
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df
import merge_pdbs


"""
Create a pdb file sample corresponding to the best scoring pdb files across all log files. This pdb file sample will be used to create a decoy dataset.
------
Params
log_files: list of Path objects to log files.
N: number of pdb files in the sample.
------
Returns
sample_df: pandas DataFrame with indexes 0-N and columns ["pdb_file"].
"""
def get_best_scoring_sample_df(
        log_files,
        N
):
    log_file_groups = list()
    log_file_groups.append(log_files)

    # Use the stat_df to get the best scoring pdb files.
    stat_df = get_stat_df.get_stat_df(
        log_file_groups=log_file_groups,
        fields=["r_free"],
        stats=["min"],
        N=N,
        offset=10,
        equil=100,
        test=False
    )

    sample_df = pd.DataFrame(index=range(N), columns=["pdb_file"])
    for i in range(N):
        log_file = stat_df.iloc[0]["r_free_min_{}_log".format(i)]
        id = stat_df.iloc[0]["r_free_min_{}_id".format(i)]

        output_dir = log_file.parents[0]
        # The index corresponds to the log file step which is 10 times the pdb file id.
        pdb_file = Path(output_dir, "pdbs", "{}.pdb".format(id//10))

        sample_df.loc[i]["pdb_file"] = pdb_file

    return sample_df


"""
Create a pdb file sample corresponding to a random sample of pdb files across all log files. This pdb file sample will be used to create a decoy dataset.
------
Params
log_files: list of Path objects to log files.
N: number of pdb files in the sample.
------
Returns
sample_df: pandas DataFrame with indexes 0-N and columns ["pdb_file"].
"""
def read_log_file(
    log_file
):
    # Log file may be empty.
    try:
        log_df = pd.read_csv(log_file, index_col=0)
    except pd.errors.EmptyDataError as e:
        return e

    return log_file, log_df


def get_random_sample_df(
        log_files,
        N
):
    # Load all log files into a single DataFrame and get a random subset.
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

        print(log_file)
        n_pdbs = (len(log_df)-100)//10
        if n_pdbs > N:
            pdb_log_df = log_df[100::10].sample(N).copy()
        else:
            pdb_log_df = log_df[100::10].sample(n_pdbs).copy()

        output_dir = log_file.parents[0]
        pdb_files = [Path(output_dir, "pdbs", "{}.pdb".format(step//10)) for step in pdb_log_df["step"]]
        pdb_log_df["pdb_file"] = pdb_files

        # To speed up the process, only include a random subset of the pdb files.
        log_dfs.append(pdb_log_df)

    merge_log_df = pd.concat(log_dfs)
    print(len(merge_log_df))
    sample_log_df = merge_log_df.sample(n=N)

    sample_df = pd.DataFrame(index=range(N), columns=["pdb_file"])
    for i in range(N):
        sample_df.loc[i]["pdb_file"] = sample_log_df.iloc[i]["pdb_file"]

    return sample_df


"""
Get a sample of pdb files from a set of log files.
------
Params
log_files: list of Path objects to log files.
N: number of pdb files in the sample.
best: if True, get the best scoring pdb files across all log files. If False, get a random sample of pdb files.
------
Returns
sample_df: pandas DataFrame with indexes 0-N and columns ["pdb_file"].
"""
def get_pdb_files(
        log_files,
        N,
        best
):
    if best:
        sample_df = get_best_scoring_sample_df(
            log_files=log_files,
            N=N
        )
        pdb_files = list(sample_df["pdb_file"])
    else:
        sample_df = get_random_sample_df(
            log_files=log_files,
            N=N
        )
        pdb_files = list(sample_df["pdb_file"])

    return pdb_files


"""
Create a decoy dataset from a sample of pdb files. The pdb file writing is multiprocessed to speed up the process.
------
Params
sample_df: pandas DataFrame with indexes 0-n_decoys*n_struct and columns ["pdb_file"].
decoy_pdb_dir: Path object to directory where decoy pdb files will be written.
n_decoys: number of decoys to create.
n_struct: number of structures in each decoy.
------
Returns
decoy_df: pandas DataFrame corresponding to all of the multi-state pdb files constructed with indexes 0-n_decoys and columns [1-n_struct].
"""
def create_decoy_dataset_from_pdb_files(
        pdb_files,
        decoy_pdb_dir,
        n_decoys,
        n_struct
):
    # There need to be enough pdb files in the sample df to populate the decoy df.
    if len(pdb_files) < n_struct*n_decoys:
        raise RuntimeError("The sample df ({}) must have enough entries to populate the decoy df ({})".format(len(sample_df), n_struct*n_decoys))

    decoy_df = pd.DataFrame(index=list(range(n_decoys)), columns=range(n_struct))

    pool_params = list()
    for i in range(n_decoys):
        out_file = Path(decoy_pdb_dir, "{}.pdb".format(i))
        decoy_pdb_files = pdb_files[i*n_struct:(i+1)*n_struct]
        for j in range(n_struct):
            decoy_df.loc[i,j]=decoy_pdb_files[j]

        params_dict = dict()
        params_dict["merge_pdb_file"] = out_file
        params_dict["pdb_files"] = decoy_pdb_files
        params_dict["occs"] = [1]*n_struct
        params_dict["n"] = -1
        params_dict["order"] = False
        params_dict["id"] = i

        pool_params.append(params_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        merge_pdbs.write_merge_pdf_file_pool,
        pool_params
    )

    for id in pool_results:
        print(id)

    return decoy_df


if __name__ == "__main__":
    job_lookup = list()
    job_lookup.append(("47_1", "3538026"))
    job_lookup.append(("48_2", "3538027"))
    job_lookup.append(("49_4", "3538028"))
    job_lookup.append(("50_8", "3538029"))
    job_lookup.append(("51_16", "3538030"))
    # job_lookup.append(("52_32", "3538031"))

    for job_pair in job_lookup:
        print(job_pair)
        job_name, job_num = job_pair

        n_struct = 1
        n_decoys = 1000
        decoy_name = "rand_1000"

        decoy_pdb_dir = Path("/wynton/group/sali/mhancock/xray/decoys/data/3ca7/{}/{}".format(job_name, decoy_name))
        decoy_pdb_dir.mkdir(exist_ok=True, parents=True)

        decoy_meta_dir = Path(Path.home(), "xray/decoys/data/3ca7/{}".format(job_name))
        decoy_meta_dir.mkdir(exist_ok=True, parents=True)
        decoy_meta_file = Path(decoy_meta_dir, "{}.csv".format(decoy_name))

        job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/{}/{}".format(job_name, job_num))

        # Get and save the dataframe containing all the decoy entries (each decoy entry containing 1 or more random pdb files and an equal number of weights).
        out_dirs = list()
        pdb_dirs = list()
        for out_dir in job_dir.glob("output_*"):
            if Path(out_dir, "pdbs").exists():
                pdb_dirs.append(Path(out_dir, "pdbs"))
                out_dirs.append(out_dir)

        # Request additional pdb files in order to ensure that there are enough valid pdb files to create the decoy dataset.
        pdb_files = get_pdb_files(
            log_files=list(job_dir.glob("output_*/log.csv")),
            N=int((n_decoys*n_struct)*1.1),
            best=False
        )
        # Check that the pdb files exist for all the entries in the sample df.
        valid_pdb_files = list()
        for pdb_file in pdb_files:
            if pdb_file.exists():
                valid_pdb_files.append(pdb_file)

        print(len(valid_pdb_files))

        # pdb_dirs = pdb_dirs[:10]
        decoy_df = create_decoy_dataset_from_pdb_files(
                pdb_files=valid_pdb_files,
                decoy_pdb_dir=decoy_pdb_dir,
                n_decoys=n_decoys,
                n_struct=n_struct
        )
        decoy_df.to_csv(decoy_meta_file)

