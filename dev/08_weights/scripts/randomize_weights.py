from pathlib import Path
import sys
import pandas as pd
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd
import multiprocessing
import random


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")

    pool_params = list()
    all_ws = list()

    for i in range(500):
        if i == 0:
            ws = [1]*4
        else:
            ws = [random.random() for j in range(4)]

        ws = [w / sum(ws) for w in ws]

        all_ws.append(ws)

        param_dict = dict()
        param_dict["decoy_file"] = Path(xray_dir, "dev/08_weights/data/0.pdb")
        param_dict["occs"] = ws
        param_dict["ref_file"] = Path(xray_dir, "data/pdbs/7mhk/7mhk_clean_h20.pdb")
        param_dict["cif_file"] = Path(xray_dir, "data/reflections/7mhk/7mhk.cif")
        param_dict["flags_file"] = Path(xray_dir, "data/reflections/7mhk/7mhk.cif")
        param_dict["res"] = 0
        param_dict["uc_dim"] = (114.300, 54.290, 44.970, 90.00, 102.12, 90.00)
        param_dict["sg"] = "C 1 2 1"
        param_dict["w_xray"] = None
        param_dict["score_fs"] = ["ml", "ff", "rmsd"]

        pool_params.append(param_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.map(
        score_rmsd.pool_score,
        pool_params
    )

    r_frees = list()
    for i in range(len(pool_results)):
        r_free = pool_results[i]["r_free"]
        r_frees.append(r_free)
        print(r_free)

    print(np.min(r_frees), np.argmin(r_frees), all_ws[np.argmin(r_frees)])


