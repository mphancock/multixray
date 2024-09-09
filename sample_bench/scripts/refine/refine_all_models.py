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

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import read_pdb_and_refine, read_pdb_and_refine_to_max_ff
from params import read_job_csv
sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
from score_rmsd import pool_score


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--job_csv_file")
    parser.add_argument("--max_ff", type=int)
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
    ref_files = params_dict["refs"]
    ref_w_mat = params_dict["ref_w_mat"]
    cif_names = [cif_file.stem for cif_file in cif_files]


    ## Get the resolution to evaluate the data at from the first entry of the sa string
    ## If multiple resolutions are used, this will need to be changed
    sample_sched = params_dict["sample_sched"]

    res = sample_sched[0]["res"]

    exp_dir = job_dir.parents[0]

    if not pdb_dir.exists():
        raise RuntimeError("PDB directory does not exist")

    log_file = Path(args.out_dir, "log.csv")
    log_df = pd.read_csv(log_file, index_col=0)
    log_df = log_df.loc[log_df["pdb"].notnull()]

    refined_exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/{}_ref_{}".format(exp_dir.name, max_ff))
    refined_out_dir = Path(refined_exp_dir, job_dir.name, out_dir.name)
    shutil.rmtree(refined_out_dir, ignore_errors=True)
    refined_pdb_dir = Path(refined_out_dir, "pdbs")
    refined_pdb_dir.mkdir(parents=True, exist_ok=True)

    refined_log_file = Path(refined_out_dir, "log.csv")
    refined_log_df = log_df.copy()

    refined_log_df["ff"] = np.nan
    refined_log_df["pdb"] = np.nan
    refined_log_df["time"] = np.nan
    for cond in range(J):
        cif_name = cif_names[cond]
        for col_type in ["r_free", "r_work", "xray"]:
            refined_log_df["{}_{}".format(col_type, cif_name)] = np.nan

    refined_log_df.reset_index(drop=True, inplace=True)
    print(len(refined_log_df))

    pdb_files = list(pdb_dir.glob("*.pdb"))
    t0 = time.time()

    # Last pdb file 301 doesn't get saved in the log file
    for i in range(len(pdb_files)-1):
    # for i in range(240, len(pdb_files)-1):
        ## Refine pdb_file
        pdb_file = pdb_files[i]
        print(pdb_file)

        if not pdb_file.exists():
            continue

        refined_pdb_file = Path(refined_pdb_dir, "{}.pdb".format(pdb_file.stem))
        print(refined_pdb_file)

        params_dict = dict()
        params_dict["pdb_file"] = pdb_file
        params_dict["out_pdb_file"] = refined_pdb_file
        params_dict["max_ff"] = max_ff
        params_dict["log_file"] = None

        read_pdb_and_refine_to_max_ff(params_dict)

        ## Now score refined_pdb_file against all cif files
        refined_log_df.loc[i, "pdb"] = refined_pdb_file
        for cond in range(J):
            cif_name = cif_names[cond]
            occs = [float(log_df.iloc[i]["w_{}_{}".format(state, cif_name)]) for state in range(N)]

            param_dict = dict()
            param_dict["decoy_file"] = refined_pdb_file
            param_dict["decoy_occs"] = occs
            param_dict["ref_file"] = ref_files[cond]
            param_dict["ref_occs"] = ref_w_mat[:, cond]
            param_dict["cif_file"] = cif_files[cond]
            param_dict["flags_file"] = cif_files[cond]
            param_dict["ab_file"] = None
            param_dict["adp_file"] = None
            param_dict["scale_k1"] = True
            param_dict["scale"] = True
            param_dict["res"] = res
            param_dict["score_fs"] = ["ml", "ff", "rmsd_avg"]

            score_dict = pool_score(param_dict)
            refined_log_df.loc[i, "xray_{}".format(cif_name)] = score_dict["ml"]
            refined_log_df.loc[i, "r_free_{}".format(cif_name)] = score_dict["r_free"]
            refined_log_df.loc[i, "r_work_{}".format(cif_name)] = score_dict["r_work"]
            refined_log_df.loc[i, "rmsd_{}".format(cif_name)] = score_dict["rmsd_avg"]
            refined_log_df.loc[i, "ff"] = score_dict["ff"]
            print(score_dict)

        refined_log_df.loc[i, "time"] = time.time() - t0

        if i % 10 == 0:
            refined_log_df.to_csv(refined_log_file)

    print(refined_log_df.head())

    refined_log_df.to_csv(refined_log_file)

