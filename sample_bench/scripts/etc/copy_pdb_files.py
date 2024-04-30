from pathlib import Path
import shutil
import math
import pandas as pd

import IMP
import IMP.atom


if __name__ == "__main__":
    summary_df = pd.read_csv(Path(Path.home(), "xray/sample_bench/data/7mhf/179_exp/summary.csv"))

    print(summary_df.head())

    for i in range(len(summary_df)):
        pdb_file = Path(summary_df.iloc[i]["pdb"])
        N = summary_df.iloc[i]["N"]
        J = summary_df.iloc[i]["J"]
        cif_name = Path(summary_df.iloc[i]["cif_name"])

        new_pdb_file = Path(Path.home(), "xray/sample_bench/data/7mhf/179_exp/summary/{}_N{}_J{}.pdb".format(cif_name.stem, N, J))

        print(pdb_file, new_pdb_file)
        shutil.copy(pdb_file, new_pdb_file)

