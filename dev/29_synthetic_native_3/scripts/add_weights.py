from pathlib import Path
import random
import numpy as np
import pandas as pd

import IMP
import IMP.atom

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
import weights


if __name__ == "__main__":
    # natives_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/7mhf_orig.csv")
    natives_file = Path(Path.home(), "xray/score_bench/data/7mhf/121_native_decoys/rand1000.csv")
    natives_df = pd.read_csv(natives_file, index_col=0)

    n_occs = 2
    n_state = 2
    floor = .05

    for i in range(len(natives_df)):
        occs_cur = None

        for j in range(n_occs):
            # pdb_file = natives_df.loc[i, "pdb"]
            # occs_cur = natives_df.loc[i, "weight_0_0"], natives_df.loc[i, "weight_0_1"]

            occs_new = weights.get_weights(
                floor=floor,
                n_state=n_state,
                occs_cur=occs_cur,
                sigma=.05
            )
            occs_cur = occs_new

            for k in range(n_state):
                natives_df.loc[i, "weight_{}_{}".format(j, k)] = occs_new[k]

    natives_file_tmp = Path(str(natives_file)+".tmp")
    natives_df.to_csv(natives_file_tmp)
