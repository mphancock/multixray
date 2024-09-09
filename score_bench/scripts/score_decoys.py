from pathlib import Path
import sys
import pandas as pd
import multiprocessing
import argparse

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd
sys.path.append(str(Path(Path.home(), "xray/src")))
import utility


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("--i", required=False)
    # parser.add_argument("--n_cond", required=True, type=int)
    parser.add_argument("--job", required=True)

    args = parser.parse_args()

    decoy_file = Path(Path.home(), "xray/score_bench/data/7mhf/{}/rand1000.csv".format(args.job))

    decoy_df = pd.read_csv(decoy_file, index_col=0)
    native_df = pd.read_csv(Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv"), index_col=0)

    n_state = utility.get_n_state_from_pdb_file(Path(decoy_df.loc[0, "pdb"]))
    ref_n_state = utility.get_n_state_from_pdb_file(Path(native_df.loc[0, "pdb"]))

    pool_params = list()
    for k in range(10):
        ref_pdb_file = Path(native_df.loc[k, "pdb"])

        for i in range(1000):
            for cond in range(2):
                cif_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/{}/{}.cif".format(k, cond))

                ref_occs = list()
                for state in range(ref_n_state):
                    ref_occs.append(native_df.loc[k, "w_{}_{}".format(state, cond)])

                # for j in range(1000):

                pdb_file = Path(decoy_df.loc[i, "pdb"])

                occs = list()
                for state in range(n_state):
                    occs.append(decoy_df.loc[i, "w_{}_{}".format(state, cond)])

                param_dict = dict()
                param_dict["decoy_file"] = pdb_file
                param_dict["decoy_occs"] = occs
                param_dict["ref_file"] = ref_pdb_file
                param_dict["ref_occs"] = ref_occs
                param_dict["cif_file"] = cif_file
                param_dict["flags_file"] = cif_file
                param_dict["res"] = 2
                param_dict["score_fs"] = ["ml", "rmsd_avg", "ff"]
                param_dict["adp_file"] = None
                param_dict["ab_file"] = None
                param_dict["scale"] = True
                param_dict["scale_k1"] = True

                pool_params.append(param_dict)

    # for pool_param in pool_params:
    #     print(pool_param)
    #     print(score_rmsd.pool_score(pool_param))

    print(multiprocessing.cpu_count())
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        score_rmsd.pool_score,
        pool_params
    )

    cnt = 0
    score_df = decoy_df.copy()
    score_df_file = Path(decoy_file.parents[0], "rand1000_score.csv")
    for score_dict in pool_results:
        # print(score_dict)
        pdb_file = score_dict["pdb_file"]
        k = score_dict["cif_file"].parents[0].name
        cond = score_dict["cif_file"].stem

        print(pdb_file, k, cond)

        score_id = score_df[score_df["pdb"] == str(pdb_file)].index[0]


        score_df.loc[score_id, "xray_{}_{}".format(cond, k)] = score_dict["ml"]
        score_df.loc[score_id, "rmsd_{}_{}".format(cond, k)] = score_dict["rmsd_avg"]
        score_df.loc[score_id, "r_free_{}_{}".format(cond, k)] = score_dict["r_free"]
        score_df.loc[score_id, "ff"] = score_dict["ff"]

        if cnt % 1000 == 0:
            score_df.to_csv(score_df_file)

        cnt = cnt + 1

    pool_obj.close()

    score_df.to_csv(score_df_file)

