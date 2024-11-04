import pickle
import pandas as pd
import numpy as np
import math
from pathlib import Path

from utility import get_n_state_from_pdb_file


def write_params_txt(
        param_dict,
        param_file
):
    param_f = open(param_file, "a")
    for key in param_dict.keys():
        print("{:<15}{}\n".format(key, param_dict[key]))
        param_f.write("{:<15}{}\n".format(key, param_dict[key]))
    param_f.write("\n\n")
    param_f.close()


def write_params_csv(
        param_dict,
        param_file
):
    param_df = pd.DataFrame()
    N = param_dict["N"]
    J = param_dict["J"]

    param_df.loc[0,"N"] = N
    param_df.loc[0,"J"] = J

    for cond in range(J):
        param_df.loc[0,"cif_{}".format(cond)] = param_dict["cifs"][cond]
        param_df.loc[0,"ref_{}".format(cond)] = param_dict["refs"][cond]

    ref_w_mat = param_dict["ref_w_mat"]
    for state in range(ref_w_mat.shape[0]):
        for cond in range(ref_w_mat.shape[1]):
            param_df.loc[0,"w_{}_{}".format(state, cond)] = ref_w_mat[state, cond]

    param_df.loc[0,"w_xray"] = param_dict["w_xray"]

    sample_sched = param_dict["sample_sched"]
    sample_sched_str = get_string_from_sample_sched(sample_sched)
    param_df.loc[0,"sample_sched_str"] = sample_sched_str

    if param_dict["start_pdb_file"]:
        param_df.loc[0,"start_pdb_file"] = param_dict["start_pdb_file"]

    param_df.to_csv(param_file)


def read_job_csv(
    job_csv_file,
    job_id
):
    job_df = pd.read_csv(job_csv_file, index_col=0)
    param_dict = dict()

    N = int(job_df.loc[job_id]["N"])
    ref_N = int(job_df.loc[job_id]["ref_N"])
    J = int(job_df.loc[job_id]["J"])

    param_dict["N"] = N
    param_dict["J"] = J

    cif_files, ref_files = list(), list()
    for cond in range(J):
        print(job_df.loc[job_id, "cif_0"], type(job_df.loc[job_id, "cif_0"]))

        if type(job_df.loc[job_id, "cif_{}".format(cond)]) == str:
            cif_files.append(Path(job_df.loc[job_id]["cif_{}".format(cond)]))
        else:
            cif_files.append(None)

    param_dict["init_pdbs"] = pdb_str_to_msmc_input(job_df.loc[job_id]["init_pdbs"], N)
    param_dict["ref_pdbs"] = pdb_str_to_msmc_input(job_df.loc[job_id]["ref_pdbs"], ref_N)

    param_dict["init_w_mat"] = build_weights_matrix(job_df, job_id, "init", N, [i for i in range(J)])
    param_dict["ref_w_mat"] = build_weights_matrix(job_df, job_id, "ref", ref_N, [i for i in range(J)])


    param_dict["cifs"] = cif_files

    for col in ["w_xray", "xray_freq", "weight_thermo", "vel_thermo", "sample_sched_str", "refine"]:
        param_dict[col] = job_df.loc[job_id][col]

    return param_dict


## msmc_m takes a list of pdb files and will read all states from all pdbs
def pdb_str_to_msmc_input(
    pdb_str,
    N
):
    ## 3 possibilities:
    ## 1) single N state pdb file
    ## 2) 1 state pdb file
    ## 3) N 1 state pdb file
    pdbs = [Path(pdb) for pdb in pdb_str.split(",")]
    pdb_n_state = get_n_state_from_pdb_file(pdbs[0])

    msmc_inputs = list()
    if pdb_n_state == N:
        msmc_inputs.append(pdbs[0])
    elif pdb_n_state == 1 and N > 1:
        for i in range(N):
            msmc_inputs.append(pdbs[0])
    else:
        for i in range(N):
            msmc_inputs.append(pdbs[i])

    return msmc_inputs


def build_weights_matrix(
    df,
    row,
    prefix,
    N,
    col_ids
):
    w_mat = np.ndarray([N, len(col_ids)])

    for state in range(N):
        for cond in range(len(col_ids)):
        # for col_id in col_ids:
            col_id = col_ids[cond]

            col = "{}_{}_{}".format(prefix, state, col_id)
            if col not in df.columns:
                raise RuntimeError("{} not in df columns.".format(col))

            w_mat[state, cond] = df.loc[row, col]

    return w_mat