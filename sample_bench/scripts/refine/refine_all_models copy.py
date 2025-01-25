from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import argparse
import time
import shutil

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import read_pdb_and_refine, read_pdb_and_refine_to_max_ff, refine_posterior
from params import read_job_csv, build_weights_matrix
sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
from score import pool_score


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--job_csv_file")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    pdb_dir = Path(out_dir, "pdbs")

    if not pdb_dir.exists():
        raise RuntimeError("{} does not exist".format(pdb_dir))

    job_dir = out_dir.parents[0]
    job_name = job_dir.name

    input_csv = Path(args.job_csv_file)
    input_id = int(job_dir.name)

    params_dict = read_job_csv(input_csv, input_id)
    N = params_dict["N"]
    J = params_dict["J"]
    cif_files = params_dict["cifs"]
    ref_files = params_dict["ref_pdbs"]
    ref_w_mat = params_dict["ref_w_mat"]
    cif_names = [cif_file.stem for cif_file in cif_files]

    print(cif_files)

    exp_dir = job_dir.parents[0]

    if not pdb_dir.exists():
        raise RuntimeError("PDB directory does not exist")

    log_file = Path(args.out_dir, "log.csv")
    log_df = pd.read_csv(log_file, index_col=0)
    log_df = log_df.loc[log_df["pdb"].notnull()]

    refined_exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/{}_phenix_ref".format(exp_dir.name))
    refined_out_dir = Path(refined_exp_dir, job_dir.name, out_dir.name)
    shutil.rmtree(refined_out_dir, ignore_errors=True)
    refined_pdb_dir = Path(refined_out_dir, "pdbs")
    refined_pdb_dir.mkdir(parents=True, exist_ok=True)

    refined_log_file = Path(refined_out_dir, "log.csv")
    refined_log_df = log_df.copy()

    refined_log_df["ff"] = np.nan
    # refined_log_df["pdb"] = np.nan
    refined_log_df["time"] = np.nan
    for cond in range(J):
        cif_name = cif_names[cond]
        for col_type in ["r_free", "r_work", "xray"]:
            refined_log_df["{}_{}".format(col_type, cif_name)] = np.nan

    ## drop any rows that have pdb as nan
    refined_log_df.reset_index(drop=True, inplace=True)

    ## pick 10 random rows to refine
    refined_log_df = refined_log_df.sample(n=10)
    refined_log_df.reset_index(drop=True, inplace=True)

    pdb_files = list(refined_log_df["pdb"])
    for i in range(len(pdb_files)):
        ## convert the pdb_file to phenix format
        #

