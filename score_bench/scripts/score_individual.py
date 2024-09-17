from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/src")))
from score import pool_score


if __name__ == "__main__":
    param_dict = dict()
    # param_dict["decoy_file"] = Path("/wynton/home/sali/mhancock/xray/dev/39_bench_ensemble/data/pdbs/7mhl.pdb")
    # param_dict["decoy_occs"] = [1/54]*54
    # param_dict["decoy_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/219_no_k1_scale/5/output_6/pdbs/500.pdb")
    param_dict["decoy_file"] = Path(Path.home(), "xray/tmp/6_500_align.pdb")
    param_dict["decoy_occs"] = [1]
    param_dict["ref_file"] = Path("/wynton/home/sali/mhancock/xray/dev/39_bench_ensemble/data/pdbs/7mhl.pdb")
    param_dict["ref_occs"] = [1/54]*54
    param_dict["cif_file"] = Path(Path.home(), "xray/dev/39_bench_ensemble/data/cifs/7mhm.cif")
    # param_dict["cif_file"] = Path(Path.home(), "xray/data/cifs/7mhf/7mhf.cif")

    param_dict["flags_file"] = param_dict["cif_file"]
    param_dict["ab_file"] = None
    param_dict["adp_file"] = None
    param_dict["scale_k1"] = True
    param_dict["scale"] = True
    param_dict["res"] = 0
    param_dict["score_fs"] = ["ml", "ff"]

    score_dict = pool_score(
        params=param_dict
    )
    print(score_dict)