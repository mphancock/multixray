from pathlib import Path
import pandas as pd
import sys

import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    log_file = Path(Path(Path.home(), "xray/tmp/log.csv"))
    ref_file = Path(Path(Path.home(), "xray/dev/17_synthetic_native/data/pdbs/2_state_uni/0.pdb"))

    n_state = 2

    log_df = pd.read_csv(log_file)

    field = "rmsd"

    pdb_log_df = log_df.dropna(subset=['pdb'])

    for i in pdb_log_df.index:
        pdb_params = dict()
        pdb_params["decoy_file"] = log_df.loc[i, "pdb"]
        pdb_params["occs"] = [log_df.loc[i, "occ_{}".format(j)] for j in range(n_state)]
        pdb_params["ref_file"] = ref_file
        pdb_params["cif_file"] = None
        pdb_params["res"] = None
        pdb_params["score_fs"] = ["rmsd_ord", "rmsd_avg"]

        score_df = score_rmsd.pool_score(pdb_params)

        print(score_df)



