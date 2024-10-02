from pathlib import Path
import sys
import pandas as pd

sys.path.append("../../src")
from score import pool_score


if __name__ == "__main__":
    param_dict = dict()
    param_dict["decoy_file"] = Path(Path.home(), "Documents/xray/dev/42_traj_analysis/data/242/pdbs/11/1.pdb")
    param_dict["decoy_occs"] = [0.5,0.5]
    param_dict["ref_file"] = Path(Path.home(), "Documents/xray/data/pdbs/3k0m/3k0m.pdb")
    param_dict["ref_occs"] = [1]
    param_dict["cif_file"] = Path(Path.home(), "Documents/xray/data/cifs/3k0m/3k0m.cif")
    param_dict["flags_file"] = param_dict["cif_file"]
    param_dict["ab_file"] = None
    param_dict["adp_type"] = "aniso"
    param_dict["adp_file"] = Path(Path.home(), "Documents/xray/data/pdbs/3k0m/3k0m.pdb")
    # param_dict["adp_type"] = "iso_avg"
    # param_dict["adp_file"] = None
    param_dict["scale_k1"] = True
    param_dict["scale"] = True
    param_dict["res"] = 0
    param_dict["score_fs"] = ["ml", "ff"]

    score_dict = pool_score(
        params=param_dict
    )
    print(score_dict)