from pathlib import Path
import sys
import pandas as pd

sys.path.append("../../src")
from score import pool_score


if __name__ == "__main__":
    param_dict = dict()
    param_dict["decoy_file"] = Path("/wynton/home/sali/mhancock/xray/dev/45_synthetic_native_4/data/pdbs/native_1.pdb")
    param_dict["decoy_occs"] = [1]
    param_dict["ref_file"] = Path("/wynton/home/sali/mhancock/xray/dev/45_synthetic_native_4/data/pdbs/native_1.pdb")
    param_dict["ref_occs"] = [1]
    param_dict["cif_file"] = Path("/wynton/home/sali/mhancock/xray/dev/45_synthetic_native_4/data/cifs/native_1.cif")
    param_dict["flags_file"] = param_dict["cif_file"]
    param_dict["ab_file"] = None
    param_dict["scale_k1"] = True
    param_dict["scale"] = True
    param_dict["res"] = 0
    param_dict["score_fs"] = ["ml", "ff"]

    score_dict = pool_score(
        params=param_dict
    )
    print(score_dict)