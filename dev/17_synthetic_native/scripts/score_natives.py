from pathlib import Path
import sys
import pandas as pd
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    rmsd_df = pd.DataFrame(index=list(range(10)))
    decoy_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs/2_state_1")
    cif_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/2_state_1_noise")
    scores_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/{}.csv".format(cif_dir.name))

    pool_params = list()
    score_fs = ["ml", "ff","rmsd_avg"]
    for i in range(10):
    # for i in [0,5]:
        pdb_file = Path(decoy_dir, "{}.pdb".format(i))
        cif_file = Path(cif_dir, "{}.cif".format(i))

        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["occs"] = None
        param_dict["ref_file"] = pdb_file
        param_dict["cif_file"] = cif_file
        param_dict["flags_file"] = cif_file
        param_dict["res"] = 2
        param_dict["score_fs"] = score_fs
        param_dict["adp_file"] = None
        param_dict["scale_k1"] = False

        pool_params.append(param_dict)

    # for pool_param in pool_params:
    #     print(pool_param)
    #     score_rmsd.pool_score(pool_param)

    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    # pool_results = pool_obj.map(
    #     score_rmsd.pool_score,
    #     pool_params
    # )

    # i = 0
    # for score_dict in pool_results:
    #     print(score_dict)
    #     for score_f in score_dict.keys():
    #         rmsd_df.loc[i, score_f] = score_dict[score_f]

    #     i = i+1

    # rmsd_df.to_csv(scores_file)

    # pool_obj.close()