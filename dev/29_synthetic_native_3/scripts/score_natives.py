from pathlib import Path
import sys
import pandas as pd
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    n_cond = 2
    n_state = 2

    native_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv")
    native_df = pd.read_csv(native_file, index_col=0)

    # for cond in list(range(n_cond)):
    pool_params = list()

    for i in range(len(native_df)):
        pdb_file = Path(native_df.loc[i, "pdb"])
        cif_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/{}".format(i))

        for cond in range(n_cond):
            occs = list()
            for state in list(range(n_state)):
                occs.append(native_df.loc[i, "w_{}_{}".format(state, cond)])

            cif_file = Path(cif_dir, "{}.cif".format(cond))

            param_dict = dict()
            param_dict["decoy_file"] = pdb_file
            param_dict["decoy_occs"] = occs
            param_dict["ref_file"] = pdb_file
            param_dict["ref_occs"] = occs
            param_dict["cif_file"] = cif_file
            param_dict["flags_file"] = cif_file
            param_dict["res"] = 2
            param_dict["score_fs"] = ["ml", "rmsd_avg", "ff"]
            param_dict["adp_file"] = None
            param_dict["ab_file"] = None
            param_dict["scale"] = True
            param_dict["scale_k1"] = True

            print(param_dict)

            pool_params.append(param_dict)

    #     # for pool_param in pool_params:
    #     #     print(pool_param)
    #     #     score_rmsd.pool_score(pool_param)

    print(multiprocessing.cpu_count())
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.map(
        score_rmsd.pool_score,
        pool_params
    )

    i = 0
    for score_dict in pool_results:
        print(score_dict)
        # native_df.loc[i, "xray_{}".format(cond)] = score_dict["ml"]
        # native_df.loc[i, "r_free_{}".format(cond)] = score_dict["r_free"]
        # native_df.loc[i, "rmsd_avg_{}".format(cond)] = score_dict["rmsd_avg"]
        # native_df.loc[i, "ff"] = score_dict["ff"]

        # i = i+1

    pool_obj.close()

    # # native_df.to_csv(str(native_file)+".tmp")
