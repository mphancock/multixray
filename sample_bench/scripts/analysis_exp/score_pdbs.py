from pathlib import Path
import pandas as pd
import sys
import multiprocessing
import argparse
import time
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
from score_rmsd import pool_score

sys.path.append(str(Path(Path.home(), "xray/src")))
from utility import get_n_state_from_pdb_file



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_file")
    args = parser.parse_args()

    sample_file = Path(args.sample_file)
    out_file = sample_file

    sample_df = pd.read_csv(sample_file, index_col=0)

    params = list()

    print(len(sample_df))

    contains_nan = sample_df['r_free'].isna().any()
    if not contains_nan:
        print("Already processed")
        exit()

    for i in range(len(sample_df)):
        pdb_file = Path(sample_df.iloc[i]["pdb"])
        if not pdb_file.exists():
            print("skipping {}".format(pdb_file))
            continue

        n_state = sample_df.iloc[i]["N"]
        cif_name = sample_df.iloc[i]["cif_name"]
        cif_file = Path(Path.home(), "xray/data/cifs/7mhf/{}.cif".format(cif_name))

        w_columns = ["w_{}".format(j) for j in range(n_state)]
        occs = [float(occ) for occ in sample_df.iloc[i][w_columns]]
        # print(occs)

        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["decoy_occs"] = occs
        param_dict["ref_file"] = Path(Path.home(), "xray/data/pdbs/7mhf/{}.pdb".format(cif_name))
        param_dict["ref_occs"] = [1]
        param_dict["cif_file"] = cif_file
        param_dict["flags_file"] = param_dict["cif_file"]
        param_dict["ab_file"] = None
        param_dict["adp_file"] = None
        param_dict["scale_k1"] = True
        param_dict["scale"] = True
        param_dict["res"] = 0
        param_dict["score_fs"] = ["ml", "ff"]
        params.append(param_dict)

    t0 = time.time()
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(pool_score, params)
    for score_dict in pool_results:
        if not score_dict:
            continue

        print(score_dict)
        index = int(Path(score_dict["pdb_file"]).stem)

        sample_df.loc[index, "r_free"] = score_dict["r_free"]
        sample_df.loc[index, "ff"] = score_dict["ff"]

    print(time.time() - t0)

    print(sample_df.head())
    sample_df.to_csv(out_file)

