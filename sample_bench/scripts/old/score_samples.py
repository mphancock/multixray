from pathlib import Path
import sys
import multiprocessing
import pandas as pd
import os
import numpy as np

import IMP
import IMP.atom
import IMP.core

sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts")))
import rmsf
sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


def pool_get_n_empty(
        params_dict
):
    job_dir = params_dict["job_dir"]
    equil = params_dict["equil"]

    n_empty = get_n_empty(
        job_dir=job_dir,
        equil=equil
    )

    return job_dir, n_empty

def get_n_empty(
        job_dir,
        equil
):
    n_empty = 0
    for out_dir in job_dir.glob("output_*"):
        pdb_dir = Path(out_dir, "pdbs")
        is_empty = False
        if len(os.listdir(pdb_dir)) < equil:
            is_empty = True
        # is_empty = dir_empty(dir_path=pdb_dir)
        if is_empty:
            n_empty = n_empty+1

    return n_empty


if __name__ == "__main__":
    glob_out_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/dev/01_w_xray/out")

    ref_pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    ref_sel = IMP.atom.AllPDBSelector()
    ref_m = IMP.Model()
    ref_h = IMP.atom.read_pdb(str(ref_pdb_file), ref_m, ref_sel)

    analysis_df = pd.DataFrame(
        columns=["job_id", "rmsd_mu", "rmsd_sig", "rmsd_min", "n_empty"],
        index=list(range(140))
    )

    n = 1000
    equil = 100

    pool_params_n_empty = list()
    for id in range(140):
        params_dict = dict()
        job_dir = Path(glob_out_dir, str(id))
        print(job_dir)

        params_dict = dict()
        params_dict["job_id"] = id
        params_dict["job_dir"] = job_dir
        params_dict["equil"] = equil
        pool_params_n_empty.append(params_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results_n_empty = pool_obj.imap(
        pool_get_n_empty,
        pool_params_n_empty
    )

    pool_params_rand_files = list()
    for job_dir, n_empty in pool_results_n_empty:
        print(job_dir, n_empty)
        job_id = int(job_dir.name)
        analysis_df.iloc[job_id, analysis_df.columns.get_loc("job_id")] = job_id
        analysis_df.iloc[job_id, analysis_df.columns.get_loc("n_empty")] = job_id

        if n_empty < 15:
            params_dict = dict()
            params_dict["job_dir"] = Path(glob_out_dir, str(job_id))
            params_dict["n"] = n
            params_dict["equil"] = 100

            pool_params_rand_files.append(params_dict)

    pool_results_rand_files = pool_obj.imap(
        rmsf.pool_get_n_random_files_from_job_dir,
        pool_params_rand_files
    )

    pool_params_score = list()
    for job_dir, sample_pdb_files in pool_results_rand_files:
        print(job_dir, sample_pdb_files)

        for pdb_file in sample_pdb_files:
            params_dict = dict()
            params_dict["decoy_file"] = pdb_file
            params_dict["ref_file"] = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
            params_dict["cif_file"] = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")
            params_dict["res"] = 4
            params_dict["uc_dim"] = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
            params_dict["sg"] = "C 1 2 1"
            params_dict["w_xray"] = 1
            params_dict["score_fs"] = ["ml", "rmsd", "rmsd_align"]

            pool_params_score.append(params_dict)

    # Second multiprocessing to get the analysis on the pdb file samples.
    pool_results_score = pool_obj.imap(
        score_rmsd.pool_score,
        pool_params_score
    )

    all_job_dict = dict()
    all_job_dict["file"] = list()
    all_job_dict["job_id"] = list()
    all_job_dict["xray_score"] = list()
    all_job_dict["rmsd"] = list()
    all_job_dict["rmsd_align"] = list()
    for scores_dict in pool_results_score:
        pdb_file = scores_dict["file"]
        xray_score = scores_dict["ml"]
        rmsd = scores_dict["rmsd"]
        rmsd_align = scores_dict["rmsd_align"]

        job_id = int(Path(pdb_file).parents[2].name)

        all_job_dict["file"].append(pdb_file)
        all_job_dict["job_id"].append(job_id)
        all_job_dict["xray_score"].append(xray_score)
        all_job_dict["rmsd"].append(rmsd)
        all_job_dict["rmsd_align"].append(rmsd_align)

    all_job_df = pd.DataFrame(all_job_dict)
    for i in range(140):
        job_df = all_job_df.loc[all_job_df["job_id"] == i]
        if len(job_df) > 0:
            job_df.to_csv(Path(Path.home(), "xray/sample_bench/dev/01_w_xray/data/sample_scores/{}.csv".format(i)))




