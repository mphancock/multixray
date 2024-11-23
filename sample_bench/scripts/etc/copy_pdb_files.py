from pathlib import Path
import shutil
import math
import pandas as pd


if __name__ == "__main__":
    job_dir = Path(Path.home(), "xray/sample_bench/data/analysis/267_full_ref")
    pdb_df_file = Path(job_dir, "summary.csv")
    pdb_df = pd.read_csv(pdb_df_file)
    print(pdb_df.head())

    pdb_dir = Path(job_dir, pdb_df_file.stem)
    pdb_dir.mkdir(exist_ok=True)

    for i in range(len(pdb_df)):
        pdb_file = Path(pdb_df.iloc[i]["pdb"])
        new_pdb_file = Path(pdb_dir, "{}.pdb".format(i))

        print(pdb_file, new_pdb_file)
        shutil.copy(pdb_file, new_pdb_file)

