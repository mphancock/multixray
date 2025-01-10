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
    parser.add_argument("--decoy_file", required=True)
    parser.add_argument("--N", required=True, type=int)

    args = parser.parse_args()

    decoy_file = Path(args.decoy_file)
    decoy_df = pd.read_csv(decoy_file, index_col=0)

    N = args.N

    pool_params = list()

    ref_pdb_file = Path(Path.home(), "xray/dev/45_synthetic_native_4/data/pdbs/native.pdb")
    ref_w_mat = np.array([[0.9],[0.1]])

    score_fs = ["w_err"]
    w_names = ["native_0"]

    for i in range(len(decoy_df)):
        ## build the w_mat
        w_mat = build_weights_matrix(decoy_df, i, "w", N, w_names)

        if "," in decoy_df.loc[i, "pdb"]:
            pdb_files = [Path(pdb) for pdb in decoy_df.loc[i, "pdb"].split(",")]
        else:
            pdb_files = [Path(decoy_df.loc[i, "pdb"])]

        param_dict = dict()
        param_dict["decoy_files"] = pdb_files
        param_dict["decoy_w_mat"] = w_mat
        param_dict["ref_file"] = ref_pdb_file
        param_dict["ref_w_mat"] = ref_w_mat
        param_dict["score_fs"] = score_fs

        pool_params.append(param_dict)

    # for pool_param in pool_params:
    #     print(pool_param)
    # print(len(pool_params))

    # pool_results = list()
    # for pool_param in pool_params:
    #     print(pool_param)
    #     pool_results.append(pool_score(pool_param))

    #     break

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
        score_df.loc[i, "w_err"] = sum(score_dict["w_err"])

        i += 1

    score_df.to_csv(score_df_file)

