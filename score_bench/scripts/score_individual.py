from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")

    param_dict = dict()

    param_dict["decoy_file"] = Path("/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb")
    param_dict["occs"] = None
    param_dict["ref_file"] = Path("/wynton/home/sali/mhancock/xray/dev/17_synthetic_native/data/pdbs/1_state_ref/0.pdb")
    param_dict["cif_file"] = Path("/wynton/home/sali/mhancock/xray/data/reflections/3ca7/3ca7_refine.cif")
    param_dict["flags_file"] = Path("/wynton/home/sali/mhancock/xray/data/reflections/3ca7_refine.cif")
    param_dict["res"] = 0
    param_dict["score_fs"] = ["ml", "rmsd_avg","rmsd_ord"]

    score_dict = score_rmsd.pool_score(
        params=param_dict
    )
    print(score_dict)