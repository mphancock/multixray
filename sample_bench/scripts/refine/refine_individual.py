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
    N = 4
    J = 2
    cif_files = [Path("/wynton/home/sali/mhancock/xray/dev/38_standard_flags/data/7mhj.cif"), Path("/wynton/home/sali/mhancock/xray/dev/38_standard_flags/data/7mhk.cif")]
    ref_files = [Path("/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_refine_15.pdb")]
    ref_w_mat = np.array([1.0])
    cif_names = [cif_file.stem for cif_file in cif_files]
    print(cif_files)

    pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/284_2_cond_4_state/14/output_0/pdbs/363.pdb")
    refined_pdb_file = Path(Path.home(), "xray/tmp/ref_out.pdb")
    print(pdb_file, refined_pdb_file)

    params_dict = dict()
    params_dict["pdb_file"] = pdb_file
    params_dict["out_pdb_file"] = refined_pdb_file

    w_mat = np.array([[0.1835615230208595, 0.3552770352176667],[0.2934976053478322,0.0332077987252963],[0.1920311647301441,0.1573457315428919],[0.3309097069011641,0.454169434514145]])

    score_dict = refine_posterior(
        pdb_file=pdb_file,
        cif_files=cif_files,
        ref_pdb_file=refined_pdb_file,
        w_mat=w_mat,
        n_step=50
    )

    #     for cif in cif_names:
    #         refined_log_df.loc[i, "xray_{}".format(cif)] = score_dict["xray_{}".format(cif)]
    #         refined_log_df.loc[i, "r_free_{}".format(cif)] = score_dict["r_free_{}".format(cif)]
    #         refined_log_df.loc[i, "r_work_{}".format(cif)] = score_dict["r_work_{}".format(cif)]

    #     refined_log_df.loc[i, "ff"] = score_dict["ff"]
    #     refined_log_df.loc[i, "pdb"] = refined_pdb_file

    # # refined_log_df.to_csv(refined_log_file)

    # for i in range(len(refined_log_df)):
    #     refined_pdb_file = Path(refined_log_df.loc[i, "pdb"])
    #     ## Now score refined_pdb_file against all cif files
    #     refined_log_df.loc[i, "pdb"] = refined_pdb_file

    #     cif_name = cif_names[cond]
    #     # occs = [float(log_df.iloc[i]["w_{}_{}".format(state, cif_name)]) for state in range(N)]
    #     decoy_w_mat = build_weights_matrix(refined_log_df, i, "w", N, cif_names)

    #     param_dict = dict()
    #     param_dict["decoy_files"] = [refined_pdb_file]
    #     param_dict["decoy_w_mat"] = decoy_w_mat

    #     ## only 1 ref file for synthetic benchmark
    #     param_dict["ref_file"] = ref_files[0]
    #     param_dict["ref_w_mat"] = ref_w_mat
    #     param_dict["score_fs"] = ["rmsd"]
    #     for cond in range(len(cif_names)):
    #         param_dict["score_fs"].append("rmsd_{}".format(cond))

    #     # param_dict["cif_file"] = cif_files[cond]
    #     # param_dict["flags_file"] = cif_files[cond]
    #     # param_dict["scale_k1"] = True
    #     # param_dict["scale"] = True
    #     # param_dict["res"] = 0

    #     score_dict = pool_score(param_dict)
    #     # refined_log_df.loc[i, "xray_{}".format(cif_name)] = score_dict["ml"]
    #     # refined_log_df.loc[i, "r_free_{}".format(cif_name)] = score_dict["r_free"]
    #     # refined_log_df.loc[i, "r_work_{}".format(cif_name)] = score_dict["r_work"]
    #     # refined_log_df.loc[i, "rmsd_{}".format(cif_name)] = score_dict["rmsd_{}".format(cif_name)]
    #     # refined_log_df.loc[i, "ff"] = score_dict["ff"]

    #     refined_log_df.loc[i, "rmsd"] = score_dict["rmsd"]
    #     for cond in range(len(cif_names)):
    #         cif_name = cif_names[cond]
    #         refined_log_df.loc[i, "rmsd_{}".format(cif_name)] = score_dict["rmsd_{}".format(cond)]

    #     print(score_dict)

    # print(refined_log_df.head())

    # refined_log_df.to_csv(refined_log_file)

