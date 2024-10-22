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
        ref_files.append(Path(job_df.loc[job_id]["ref_{}".format(cond)]))

    ref_n_state = get_n_state_from_pdb_file(ref_files[0])
    ref_w_mat = np.ndarray([ref_n_state, J])

    # If the job csv file contains reference weights, use them. Otherwise, use uniform weights.
    rows = str(job_df.loc[job_id]["ref_weights"]).split(';')
    matrix = [list(map(float, row.split(','))) for row in rows]
    ref_w_mat = np.array(matrix)
    ref_w_mat = ref_w_mat.T
    # for cond in range(J):
        # for state in range(ref_n_state):
        #     col = "w_{}_{}".format(state, cond)
        #     if col in job_df.columns:
        #         ref_w_mat[state, cond] = job_df.loc[job_id]["w_{}_{}".format(state, cond)]
        #     else:
        #         ref_w_mat[state, cond] = 1/ref_n_state

    rows = str(job_df.loc[job_id]["init_weights"]).split(';')
    matrix = [list(map(float, row.split(','))) for row in rows]
    w_mat = np.array(matrix)
    w_mat = w_mat.T
    # w_mat = np.ndarray([N, J])
    # for cond in range(J):
    #     if job_df.loc[job_id]["init_weights"] == "uni":
    #         w_mat[:, cond] = 1/N
    #     elif job_df.loc[job_id]["init_weights"]:
    #         rows = job_df.loc[job_id]["init_weights"].split(';')
    #         matrix = [list(map(float, row.split(','))) for row in rows]
    #         w_mat = np.array(matrix)
    #     else:
    #         raise RuntimeError("No initial weights provided.")
        # else:
        #     w_mat[:, cond] = weights.get_weights(floor=0.05, n_state=N)
    param_dict["w_mat"] = w_mat

    param_dict["cifs"] = cif_files
    param_dict["refs"] = ref_files
    param_dict["ref_w_mat"] = ref_w_mat

    for col in ["w_xray", "xray_freq", "weight_thermo", "vel_thermo", "sample_sched_str", "refine"]:
        param_dict[col] = job_df.loc[job_id][col]

    start_pdb_files = job_df.loc[job_id]["start_pdb_file"].split(",")
    start_pdb_files = [Path(pdb_file) for pdb_file in start_pdb_files]
    n_state_pdb_file = get_n_state_from_pdb_file(start_pdb_files[0])

    if n_state_pdb_file == N:
        param_dict["start_pdb_file"] = [start_pdb_files[0]]
    elif n_state_pdb_file == 1 and N > 1 and len(start_pdb_files) == 1:
        param_dict["start_pdb_file"] = [Path(start_pdb_files[0])]*N
    else:
        param_dict["start_pdb_file"] = start_pdb_files

    # Optional params
    for param in ["init_weights"]:
        if param in job_df.columns:
            param_dict[param] = job_df.loc[job_id, param]
        else:
            param_dict[param] = None

    return param_dict
