from pathlib import Path
import sys
import pandas as pd
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    rmsd_df = pd.DataFrame(index=list(range(40)))

    for decoy_name in ["1_state_ref", "2_state_ref", "2_state_uni", "4_state_ref"]:
        scores_file = Path(Path.home(), "xray/dev/17_synthetic_native/data/scores/{}.csv".format(decoy_name))
        pool_params = list()

        score_fs = ["ml", "ff","rmsd_avg","rmsd_ord","rmsd_dom", "weight_delta"]

        for i in range(40):
            pdb_file = Path(Path.home(), "xray/dev/17_synthetic_native/data/pdbs/{}/{}.pdb".format(decoy_name, i))

            cif_file = Path(Path.home(), "xray/dev/17_synthetic_native/data/cifs/{}/{}.cif".format(decoy_name, i))

            param_dict = dict()
            param_dict["decoy_file"] = Path("/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb")
            param_dict["occs"] = [1]
            param_dict["ref_file"] = pdb_file
            param_dict["cif_file"] = cif_file
            param_dict["flags_file"] = cif_file
            param_dict["res"] = 0
            param_dict["score_fs"] = score_fs

            pool_params.append(param_dict)

        # for pool_param in pool_params:
        #     print(pool_param)
        #     score_rmsd.pool_score(pool_param)

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
            for score_f in score_dict.keys():
                rmsd_df.loc[i, score_f] = score_dict[score_f]

            i = i+1

        rmsd_df.to_csv(scores_file)

        pool_obj.close()