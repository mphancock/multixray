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
    parser.add_argument("--max_ff", type=int)
    parser.add_argument("--end_only", action="store_true")
    args = parser.parse_args()

    max_ff = args.max_ff

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

    refined_exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/{}_ref".format(exp_dir.name))
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

    ## if end only then only refine the last row
    if args.end_only:
        ## need to [[ ]] else it will return a series
        refined_log_df = refined_log_df.loc[[len(refined_log_df)-1]]
        refined_log_df.reset_index(drop=True, inplace=True)
    ## else pick 10 random rows to refine
    else:
        refined_log_df = refined_log_df.sample(n=10)
        refined_log_df.reset_index(drop=True, inplace=True)

    pdb_files = list(refined_log_df["pdb"])
    for i in range(len(pdb_files)):
        ## Refine pdb_file
        pdb_file = Path(pdb_files[i])

        if not pdb_file.exists():
            continue

        refined_pdb_file = Path(refined_pdb_dir, "{}.pdb".format(pdb_file.stem))
        print(pdb_file, refined_pdb_file)

        params_dict = dict()
        params_dict["pdb_file"] = pdb_file
        params_dict["out_pdb_file"] = refined_pdb_file

        w_mat = build_weights_matrix(refined_log_df, i, "w", N, cif_names)

        score_dict = refine_posterior(
            pdb_file=pdb_file,
            cif_files=cif_files,
            ref_pdb_file=refined_pdb_file,
            w_mat=w_mat,
            n_step=50
        )

        for cif in cif_names:
            refined_log_df.loc[i, "xray_{}".format(cif)] = score_dict["xray_{}".format(cif)]
            refined_log_df.loc[i, "r_free_{}".format(cif)] = score_dict["r_free_{}".format(cif)]
            refined_log_df.loc[i, "r_work_{}".format(cif)] = score_dict["r_work_{}".format(cif)]

        refined_log_df.loc[i, "ff"] = score_dict["ff"]
        refined_log_df.loc[i, "pdb"] = refined_pdb_file

    # refined_log_df.to_csv(refined_log_file)

    for i in range(len(refined_log_df)):
        refined_pdb_file = Path(refined_log_df.loc[i, "pdb"])
        ## Now score refined_pdb_file against all cif files
        refined_log_df.loc[i, "pdb"] = refined_pdb_file

        cif_name = cif_names[cond]
        # occs = [float(log_df.iloc[i]["w_{}_{}".format(state, cif_name)]) for state in range(N)]
        decoy_w_mat = build_weights_matrix(refined_log_df, i, "w", N, cif_names)

        param_dict = dict()
        param_dict["decoy_files"] = [refined_pdb_file]
        param_dict["decoy_w_mat"] = decoy_w_mat

        ## only 1 ref file for synthetic benchmark
        param_dict["ref_file"] = ref_files[0]
        param_dict["ref_w_mat"] = ref_w_mat
        param_dict["score_fs"] = ["rmsd"]
        for cond in range(len(cif_names)):
            param_dict["score_fs"].append("rmsd_{}".format(cond))

        # param_dict["cif_file"] = cif_files[cond]
        # param_dict["flags_file"] = cif_files[cond]
        # param_dict["scale_k1"] = True
        # param_dict["scale"] = True
        # param_dict["res"] = 0

        score_dict = pool_score(param_dict)
        # refined_log_df.loc[i, "xray_{}".format(cif_name)] = score_dict["ml"]
        # refined_log_df.loc[i, "r_free_{}".format(cif_name)] = score_dict["r_free"]
        # refined_log_df.loc[i, "r_work_{}".format(cif_name)] = score_dict["r_work"]
        # refined_log_df.loc[i, "rmsd_{}".format(cif_name)] = score_dict["rmsd_{}".format(cif_name)]
        # refined_log_df.loc[i, "ff"] = score_dict["ff"]

        refined_log_df.loc[i, "rmsd"] = score_dict["rmsd"]
        for cond in range(len(cif_names)):
            cif_name = cif_names[cond]
            refined_log_df.loc[i, "rmsd_{}".format(cif_name)] = score_dict["rmsd_{}".format(cond)]

        print(score_dict)

    print(refined_log_df.head())

    refined_log_df.to_csv(refined_log_file)

