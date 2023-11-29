from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    data_dir = Path(Path.home(), "xray/dev/19_synthetic_native_2/data")
    score_dir = Path(data_dir, "scores")

    pdb_dirs, cif_dirs = list(), list()
    pdb_dirs.append(Path(data_dir, "pdbs/2_state_0"))
    cif_dirs.append(Path(data_dir, "cifs/2_state_0"))
    # for i in range(2,5):
    #     pdb_dirs.append(Path(data_dir, "pdbs/4_state_{}".format(i)))
    #     cif_dirs.append(Path(data_dir, "cifs/4_state_{}".format(i)))

    n_native = 10
    n_state = 2

    for i in range(len(pdb_dirs)):
        pdb_dir = pdb_dirs[i]
        cif_dir = cif_dirs[i]
        score_file = Path(score_dir, "{}.csv".format(pdb_dir.name))
        score_df = pd.DataFrame()

        for j in range(n_native):
            native_pdb_file = Path(pdb_dir, "{}.pdb".format(j))
            native_cif_file = Path(cif_dir, "{}.cif".format(j))

            param_dict = dict()
            param_dict["decoy_file"] = native_pdb_file
            param_dict["occs"] = None
            param_dict["ref_file"] = native_pdb_file
            param_dict["cif_file"] = native_cif_file
            param_dict["flags_file"] = param_dict["cif_file"]
            param_dict["adp_file"] = None
            param_dict["res"] = 0
            param_dict["score_fs"] = ["ml", "rmsd_avg", "ff", "w"]

            score_dict = score_rmsd.pool_score(
                params=param_dict
            )
            print(score_dict)

            score_df.loc[j, "pdb"] = str(native_pdb_file)
            score_df.loc[j, "xray_0"] = score_dict["ml"]
            score_df.loc[j, "r_free_0"] = score_dict["r_free"]
            score_df.loc[j, "rmsd_avg_0"] = score_dict["rmsd_avg"]
            score_df.loc[j, "ff"] = score_dict["ff"]

            for k in range(n_state):
                score_df.loc[j, "weight_0_{}".format(k)] = score_dict["w"][k]

        score_df.to_csv(score_file)



