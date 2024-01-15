from pathlib import Path
import sys
import pandas as pd
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    pdb_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs")
    cif_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/1")
    pdb_files = list(pdb_dir.glob("*.pdb"))
    w_set = 1

    native_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/natives.csv")
    native_df = pd.read_csv(native_file, index_col=0)

    pool_params = list()
    for i in range(len(pdb_files)):
        pdb_file = pdb_files[i]
        w_0 = native_df.loc[int(pdb_file.stem), "weight_{}_0".format(w_set)]
        w_1 = native_df.loc[int(pdb_file.stem), "weight_{}_1".format(w_set)]

        cif_file = Path(cif_dir, "{}.cif".format(pdb_file.stem))

        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["occs"] = (w_0, w_1)
        param_dict["ref_file"] = pdb_file
        param_dict["cif_file"] = cif_file
        param_dict["flags_file"] = cif_file
        param_dict["res"] = 2
        param_dict["score_fs"] = ["ml", "rmsd_avg", "ff"]
        param_dict["adp_file"] = None
        param_dict["scale_k1"] = False

        pool_params.append(param_dict)

    # for pool_param in pool_params:
    #     print(pool_param)
    #     score_rmsd.pool_score(pool_param)

    print(multiprocessing.cpu_count())
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.map(
        score_rmsd.pool_score,
        pool_params
    )

    # i = 0
    for score_dict in pool_results:
        print(score_dict)

        pdb_id = int(Path(score_dict["pdb_file"]).stem)

        native_df.loc[pdb_id, "xray_{}".format(w_set)] = score_dict["ml"]
        native_df.loc[pdb_id, "rmsd_{}".format(w_set)] = score_dict["rmsd_avg"]
        native_df.loc[pdb_id, "ff"] = score_dict["ff"]

    print(native_df.head())

    #     i = i+1

    native_df.to_csv(str(native_file)+".tmp")

    pool_obj.close()