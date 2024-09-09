from pathlib import Path
import shutil
import math
import pandas as pd

import IMP
import IMP.atom


if __name__ == "__main__":
    summary_df = pd.read_csv(Path(Path.home(), "xray/sample_bench/data/7mhf/189_exp_ref_10000/sample.csv"))

    print(summary_df.head())

    for i in range(len(summary_df)):
        pdb_file = Path(summary_df.iloc[i]["pdb"])
        # N = summary_df.iloc[i]["N"]
        # J = summary_df.iloc[i]["J"]
        # cif_name = Path(summary_df.iloc[i]["cif_name"])

        new_pdb_file = Path(Path.home(), "xray/sample_bench/data/7mhf/189_exp_ref_10000/sample/{}.pdb".format(i))

        print(pdb_file, new_pdb_file)
        shutil.copy(pdb_file, new_pdb_file)

