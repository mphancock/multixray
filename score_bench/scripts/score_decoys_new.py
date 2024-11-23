from pathlib import Path
import sys
import pandas as pd
import multiprocessing
import argparse
import numpy as np

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from score import pool_score
from params import build_weights_matrix


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("--i", required=False)
    # parser.add_argument("--n_cond", required=True, type=int)
    parser.add_argument("--decoy_file", required=True)

    args = parser.parse_args()

    # decoy_file = Path(Path.home(), "xray/dev/48_decoys/data/csvs/decoys_w.csv")
    decoy_file = Path(args.decoy_file)

    decoy_df = pd.read_csv(decoy_file, index_col=0)

    N, J = 1, 1

    pool_params = list()

    cif_files = [Path(Path.home(), "xray/dev/45_synthetic_native_4/data/cifs/native_0.cif")]

    ref_pdb_file = Path(Path.home(), "xray/dev/45_synthetic_native_4/data/pdbs/native.pdb")
    ref_w_mat = np.array([[0.9],[0.1]])

    for i in range(len(decoy_df)):
        ## build the w_mat
        w_mat = build_weights_matrix(decoy_df, i, "w", N, [str(cond) for cond in range(J)])

        pdb_file = Path(decoy_df.loc[i, "pdb"])

        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["decoy_w_mat"] = w_mat
        param_dict["ref_file"] = ref_pdb_file
        param_dict["ref_w_mat"] = ref_w_mat
        param_dict["score_fs"] = ["xray_0", "rmsd_0", "ff"]

        param_dict["cif_files"] = cif_files
        param_dict["res"] = 0
        param_dict["scale"] = True
        param_dict["scale_k1"] = True

        pool_params.append(param_dict)

    # for pool_param in pool_params:
    #     print(pool_param)
    # print(len(pool_params))

    # pool_results = list()
    # for pool_param in pool_params:
    #     print(pool_param)
    #     pool_results.append(pool_score(pool_param))

    print(multiprocessing.cpu_count())
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_score,
        pool_params
    )

    i = 0
    score_df = decoy_df.copy()
    score_df_file = Path(decoy_file.parents[0], "{}_score.csv".format(decoy_file.stem))
    for score_dict in pool_results:
        # print(score_dict)
        cif_files = score_dict["cif_files"]

        pdb_file = score_dict["pdb_file"]

        for j in range(len(cif_files)):
            cif_name = cif_files[j].stem
            print(pdb_file, cif_name)

            score_df.loc[i, "xray_{}".format(cif_name)] = score_dict["xray_{}".format(j)]
            score_df.loc[i, "rmsd_{}".format(cif_name)] = score_dict["rmsd_{}".format(j)]
            score_df.loc[i, "r_free_{}".format(cif_name)] = score_dict["r_free_{}".format(j)]
            score_df.loc[i, "r_work_{}".format(cif_name)] = score_dict["r_work_{}".format(j)]
            score_df.loc[i, "ff"] = score_dict["ff"]

        i += 1

    score_df.to_csv(score_df_file)

