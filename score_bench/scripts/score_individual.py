from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")

    param_dict = dict()
    # param_dict["decoy_file"] = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7.pdb")
    # param_dict["occs"] = [1]*1
    # param_dict["ref_file"] = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7.pdb")
    # param_dict["cif_file"] = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")
    # param_dict["flags_file"] = Path(xray_dir, "data/reflections/3ca7/3ca7.cif")
    # param_dict["res"] = 0
    # param_dict["uc_dim"] = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    # param_dict["sg"] = "C 1 2 1"
    # param_dict["w_xray"] = 30000
    # param_dict["score_fs"] = ["ml", "ff", "rmsd"]

    param_dict["decoy_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/test/output_0/pdbs/1.pdb")
    param_dict["occs"] = [1]*1
    param_dict["ref_file"] = Path(xray_dir, "data/pdbs/7mhk/7mhk_clean.pdb")
    param_dict["cif_file"] = Path(xray_dir, "data/reflections/7mhk/7mhk.cif")
    param_dict["flags_file"] = Path(xray_dir, "data/reflections/7mhk/7mhk.cif")
    param_dict["res"] = 0
    param_dict["uc_dim"] = (114.300, 54.290, 44.970, 90.00, 102.12, 90.00)
    param_dict["sg"] = "C 1 2 1"
    param_dict["w_xray"] = None
    param_dict["score_fs"] = ["ml", "ff", "rmsd"]

    score_dict = score_rmsd.pool_score(
        params=param_dict
    )
    print(score_dict)