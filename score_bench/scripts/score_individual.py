from pathlib import Path
import sys
import pandas as pd
import numpy as np

sys.path.append("../../src")
from score import pool_score


if __name__ == "__main__":
    param_dict = dict()
    param_dict["decoy_files"] = [Path("/wynton/home/sali/mhancock/xray/tmp/7mhf_clean.pdb")]
    param_dict["decoy_w_mat"] = np.array([[1]])
    param_dict["ref_file"] = Path("/wynton/home/sali/mhancock/xray/dev/45_synthetic_native_4/data/pdbs/native.pdb")
    param_dict["ref_w_mat"] = np.array([[0.5], [0.5]])
    param_dict["cif_files"] = [Path("/wynton/home/sali/mhancock/xray/dev/38_standard_flags/data/7mhi.cif")]
    # param_dict["flags_file"] = param_dict["cif_fil    e"]
    param_dict["ab_file"] = None
    param_dict["scale_k1"] = True
    param_dict["scale"] = True
    param_dict["remove_outliers"] = False
    param_dict["res"] = 0
    param_dict["score_fs"] = ["xray_0", "ff"]

    score_dict = pool_score(
        params=param_dict
    )
    print(score_dict)