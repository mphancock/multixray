from pathlib import Path
import pandas as pd
import sys
import multiprocessing
import argparse
import time

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
from score_rmsd import pool_score

sys.path.append(str(Path(Path.home(), "xray/src")))
from utility import get_n_state_from_pdb_file



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_df_file")
    parser.add_argument("--skip", required=False, type=int)
    args = parser.parse_args()

    pdb_df_file = Path(args.pdb_df_file)

    # pdb_df_file = Path(Path.home(), "xray/sample_bench/data/7mhf/{}/ref.csv".format(exp_name))
    pdb_df = pd.read_csv(pdb_df_file, index_col=0)

    params = list()

    if args.skip:
        start = args.skip
    else:
        start = 0

    # for i in [14833]:
    print(len(pdb_df))
    for i in range(start, len(pdb_df)):
        pdb_file = Path(pdb_df.loc[i, "pdb"])
        if not pdb_file.exists():
            print("skipping {}".format(pdb_file))
            continue

        n_state = pdb_df.loc[i, "N"]
        cif_name = pdb_df.loc[i, "cif_name"]
        cif_file = Path(Path.home(), "xray/data/cifs/7mhf/{}.cif".format(cif_name))

        w_columns = ["w_{}".format(j) for j in range(n_state)]
        occs = [float(occ) for occ in pdb_df.loc[i, w_columns]]
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

        # print(params)

    t0 = time.time()
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(pool_score, params)
    j = 0
    for score_dict in pool_results:
        if not score_dict:
            continue

        print(score_dict)

        i = int(Path(score_dict["pdb_file"]).stem)
        pdb_df.loc[i, "r_free"] = score_dict["r_free"]
        pdb_df.loc[i, "ff"] = score_dict["ff"]

        if j % 100 == 0:
            pdb_df.to_csv(str(pdb_df_file)+".tmp")

        j = j + 1

    print(time.time() - t0)

    pdb_df.to_csv(str(pdb_df_file)+".tmp")

