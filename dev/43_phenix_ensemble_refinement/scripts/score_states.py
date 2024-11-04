from pathlib import Path
import multiprocessing

import pandas as pd

import IMP
import IMP.atom

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
# from cctbx_score import get_score
from score import pool_score
from miller_ops import get_miller_array, clean_miller_array

if __name__ == "__main__":
    pool_params = list()

    pdb_name = "3k0n"
    scores_file = Path("../data/scores/{}.csv".format(pdb_name))

    pdb_dir = Path("../data/pdbs_state/{}".format(pdb_name))
    pdb_files = list(pdb_dir.glob("*.pdb"))
    for pdb_file in pdb_files:
        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["decoy_occs"] = [1]
        param_dict["ref_file"] = Path(Path.home(), "xray/data/pdbs/3k0m/3k0m_clean.pdb")
        param_dict["ref_occs"] = [1]
        param_dict["cif_file"] = Path(Path.home(), "xray/data/cifs/3k0m/{}.cif".format(pdb_name))
        param_dict["res"] = 0
        param_dict["score_fs"] = ["ml", "rmsd_avg", "ff"]
        param_dict["adp_file"] = None
        param_dict["scale_k1"] = True
        param_dict["scale"] = True

        pool_params.append(param_dict)

    # for param in pool_params:
    #     print(pool_score(param))

    #     break

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_score,
        pool_params
    )

    # columns = ["pdb_file", "native"]
    # columns.extend(score_fs)
    # columns = score_dict.keys()
    all_scores_df = pd.DataFrame()
    i = 0
    for score_dict in pool_results:
        print(score_dict)
        for key in score_dict.keys():
            all_scores_df.loc[i, key] = score_dict[key]

        i = i+1

    all_scores_df.to_csv(scores_file)