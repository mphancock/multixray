from pathlib import Path
import pandas as pd
import shutil


if __name__ == "__main__":
    # native_df = pd.read_csv(Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv"))
    native_df = pd.read_csv(Path(Path.home(), "xray/sample_bench/data/7mhf/141_native_4_state_1_cond/score_analysis_rmsd_0.csv"))

    for i in range(10):
        # pdb_file = native_df.loc[i, "pdb"]
        pdb_file = native_df.loc[i, "pdb"]
        print(pdb_file)

        # shutil.copy(pdb_file, Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30/{}.pdb".format(i)))
        shutil.copy(pdb_file, Path(Path.home(), "xray/dev/33_best_scoring/data/4_state/{}.pdb".format(i)))

