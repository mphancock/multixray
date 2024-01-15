from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")

    param_dict = dict()

    # param_dict["decoy_file"] = Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/0.pdb")
    param_dict["decoy_file"] = Path("/wynton/home/sali/mhancock/xray/dev/32_score/data/decoy.pdb")
    param_dict["occs"] = None
    param_dict["ref_file"] = param_dict["decoy_file"]
    param_dict["cif_file"] = Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/0/0.cif")
    param_dict["flags_file"] = param_dict["cif_file"]
    # param_dict["adp_file"] = Path("/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/adps/1_state/3ca7_refine_b_factor_refine_001.pdb")
    param_dict["adp_file"] = None
    param_dict["res"] = 2
    param_dict["score_fs"] = ["ml", "rmsd_avg", "ff"]
    param_dict["scale_k1"] = False

    score_dict = score_rmsd.pool_score(
        params=param_dict
    )
    print(score_dict)