from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")

    param_dict = dict()

    # param_dict["decoy_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/87_native_4x/0/output_375/pdbs/572.pdb")
    param_dict["decoy_file"] = Path("/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/pdbs/4_state/5.pdb")

    param_dict["occs"] = None
    param_dict["ref_file"] = Path("/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/pdbs/4_state/5.pdb")
    param_dict["cif_file"] = Path("/wynton/home/sali/mhancock/xray/tmp/5.cif")
    param_dict["flags_file"] = param_dict["cif_file"]
    # param_dict["adp_file"] = Path("/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/adps/1_state/3ca7_refine_b_factor_refine_001.pdb")

    # param_dict["adp_file"] = Path("/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/adps/4_state/3ca7_refine_b_factor_refine_001.pdb")
    param_dict["adp_file"] = None
    param_dict["res"] = 0
    param_dict["score_fs"] = ["ml", "rmsd_avg"]

    score_dict = score_rmsd.pool_score(
        params=param_dict
    )
    print(score_dict)