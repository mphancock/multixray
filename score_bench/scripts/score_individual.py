from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    param_dict = dict()
    param_dict["decoy_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/test/output_0/pdbs/10.pdb")
    # param_dict["decoy_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/179_exp/N0_J0/output_8/pdbs/104.pdb")
    param_dict["decoy_occs"] = [0.7871057220298390,0.2128942779701610]
    param_dict["ref_file"] = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30/0.pdb")
    param_dict["ref_occs"] = [0.5,0.5]
    param_dict["cif_file"] = Path("/wynton/home/sali/mhancock/xray/tmp/0.cif")
    # param_dict["cif_file"] = Path(Path.home(), "xray/data/cifs/7mhf/7mhf.cif")
    param_dict["flags_file"] = param_dict["cif_file"]
    param_dict["ab_file"] = None
    param_dict["adp_file"] = None
    param_dict["scale_k1"] = True
    param_dict["scale"] = True
    param_dict["res"] = None
    param_dict["score_fs"] = ["ml", "ff"]

    score_dict = score_rmsd.pool_score(
        params=param_dict
    )
    print(score_dict)