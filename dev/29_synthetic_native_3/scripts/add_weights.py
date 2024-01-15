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
    natives_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/natives.csv")
    natives_df = pd.read_csv(natives_file, index_col=0)

    for i in range(len(natives_df)):
        pdb_file = natives_df.loc[i, "pdb"]
        ws = natives_df.loc[i, "weight_0_0"], natives_df.loc[i, "weight_0_1"]

        ws_new = weights.get_weights(
            floor=.05,
            ws_cur=ws,
            sigma=.05
        )

        print(pdb_file, ws, ws_new)

        natives_df.loc[i, "weight_1_0"] = ws_new[0]
        natives_df.loc[i, "weight_1_1"] = ws_new[1]

    natives_file_tmp = Path(str(natives_file)+".tmp")
    natives_df.to_csv(natives_file_tmp)
