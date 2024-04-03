from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")

    param_dict = dict()

    param_dict["decoy_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/169_N8/55/output_349/pdbs/277.pdb")
    param_dict["decoy_occs"] = [0.1707090564501260, 0.2154425789446470,	0.135708096184682, 0.0539801939285565, 0.1220519934501550, 0.1262063697376800, 0.0563696183735114, 0.1195320929306410]
    param_dict["ref_file"] = Path("/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi.pdb")
    param_dict["ref_occs"] = [1]
    param_dict["cif_file"] = Path("/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhi.cif")
    param_dict["flags_file"] = param_dict["cif_file"]
    param_dict["ab_file"] = None
    param_dict["adp_file"] = None
    param_dict["scale_k1"] = True
    param_dict["scale"] = True
    param_dict["res"] = 0
    param_dict["score_fs"] = ["ml", "ff"]

    score_dict = score_rmsd.pool_score(
        params=param_dict
    )
    print(score_dict)