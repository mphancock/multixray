from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")

    param_dict = dict()

    param_dict["decoy_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/121_native_decoys/7/output_620/pdbs/16.pdb")
    param_dict["decoy_w"] = [0.3747672651241047,0.6252327348758953]
    param_dict["ref_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/54_7mhf_decoys_100/4585081/output_7/pdbs/30.pdb")
    param_dict["ref_w"] = [0.297701472338339,0.702298527661661]
    param_dict["cif_file"] = Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/0/0.cif")
    param_dict["flags_file"] = param_dict["cif_file"]
    param_dict["ab_file"] = None
    param_dict["adp_file"] = None
    param_dict["scale_k1"] = True
    param_dict["scale"] = True
    param_dict["res"] = 2
    param_dict["score_fs"] = ["ml", "rmsd_avg", "ff"]

    score_dict = score_rmsd.pool_score(
        params=param_dict
    )
    print(score_dict)