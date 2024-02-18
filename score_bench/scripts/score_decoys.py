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
    parser.add_argument("--i", required=False)
    parser.add_argument("--n_cond", required=True, type=int)
    parser.add_argument("--job", required=True)

    args = parser.parse_args()

    decoy_file = Path(Path.home(), "xray/score_bench/data/7mhf/{}/rand1000.csv".format(args.job))

    decoy_df = pd.read_csv(decoy_file, index_col=0)
    native_df = pd.read_csv(Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv"), index_col=0)

    if args.i is not None:
        i_set = [int(args.i)]
    else:
        i_set = list(range(10))

    n_state = utility.get_n_state_from_pdb_file(Path(decoy_df.loc[0, "pdb"]))
    ref_n_state = utility.get_n_state_from_pdb_file(Path(native_df.loc[0, "pdb"]))

    for i in i_set:
        for cond in list(range(args.n_cond)):
            cif_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/{}/{}.cif".format(cond, i))

            ref_occs = list()
            for state in range(ref_n_state):
                ref_occs.append(native_df.loc[i, "w_{}_{}".format(state, cond)])

            pool_params = list()
            for j in range(len(decoy_df)):
            # for j in range(1000):
                pdb_file = decoy_df.loc[j, "pdb"]

                occs = list()
                for state in range(n_state):
                    occs.append(decoy_df.loc[j, "w_{}_{}".format(state, cond)])

                param_dict = dict()
                param_dict["decoy_file"] = pdb_file
                param_dict["decoy_w"] = occs
                param_dict["ref_file"] = Path(native_df.loc[i, "pdb"])
                param_dict["ref_w"] = ref_occs
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
            #     score_rmsd.pool_score(pool_param)

            print(multiprocessing.cpu_count())
            pool_obj = multiprocessing.Pool(
                multiprocessing.cpu_count()
            )

            pool_results = pool_obj.imap(
                score_rmsd.pool_score,
                pool_params
            )

            j = 0
            for score_dict in pool_results:
                print(score_dict)
                decoy_df.loc[j, "xray_{}".format(cond)] = score_dict["ml"]
                decoy_df.loc[j, "rmsd_{}".format(cond)] = score_dict["rmsd_avg"]
                decoy_df.loc[j, "r_free_{}".format(cond)] = score_dict["r_free"]
                # decoy_df.loc[j, "ff"] = score_dict["ff"]

                j = j+1

            pool_obj.close()

        decoy_df.to_csv(Path(decoy_file.parents[0], "rand1000_{}.csv".format(i)))

