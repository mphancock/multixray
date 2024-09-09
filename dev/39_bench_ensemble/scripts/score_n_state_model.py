from pathlib import Path
import sys
import multiprocessing

import pandas as pd
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/src")))
from score import pool_score
from merge_pdbs import write_merge_pdb_file



if __name__ == "__main__":
    # pdb_1_state_files = [pdb_file for pdb_file in Path("../data/pdbs/7mhl").glob("*.pdb")]
    # print(pdb_1_state_files)

    # for i in range(1,55):
    #     print(i)
    #     pdb_n_state_dir = Path("../data/pdbs/7mhl_n_state/{}".format(i))
    #     pdb_n_state_dir.mkdir(exist_ok=True)

    #     rand_n_pdb_files = np.random.choice(pdb_1_state_files, i, replace=False)

    #     write_merge_pdb_file(
    #         merge_pdb_file=Path(pdb_n_state_dir, "0.pdb"),
    #         pdb_files=rand_n_pdb_files,
    #         occs=[1]*i,
    #         n=-1,
    #         order=False,
    #         state=0
    #     )

    params = list()
    for i in range(1,55):
        pdb_file = Path("../data/pdbs/7mhl_n_state/{}/0.pdb".format(i))
        print(pdb_file)
        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["decoy_occs"] = [1/i]*i
        param_dict["ref_file"] = Path("/wynton/home/sali/mhancock/xray/dev/39_bench_ensemble/data/pdbs/7mhl.pdb")
        param_dict["ref_occs"] = [1/54]*54
        param_dict["cif_file"] = Path(Path.home(), "xray/dev/39_bench_ensemble/data/cifs/7mhl.cif")
        param_dict["flags_file"] = param_dict["cif_file"]
        param_dict["ab_file"] = None
        param_dict["adp_file"] = None
        param_dict["scale_k1"] = True
        param_dict["scale"] = True
        param_dict["res"] = 0
        param_dict["score_fs"] = ["ml", "ff"]

        params.append(param_dict)

    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_score,
        params
    )

    r_frees, xrays, pdb_ids = list(), list(), list()
    for score_dict in pool_results:
        print(score_dict)
        r_free = score_dict["r_free"]
        ml = score_dict["ml"]
        pdb_file = score_dict["pdb_file"]
        n_state = int(Path(pdb_file).parents[0].stem)

        r_frees.append(r_free)
        xrays.append(ml)
        pdb_ids.append(n_state)

    scores_df = pd.DataFrame({
        "pdb_id": pdb_ids,
        "r_free": r_frees,
        "ml": xrays
    })
    scores_df.to_csv("../data/scores/n_state.csv")

    # print(np.mean(r_frees), np.mean(xrays))


