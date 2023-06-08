from pathlib import Path
import pandas as pd
import shutil


if __name__ == "__main__":
    decoy_meta_file = Path(Path.home(), "xray/decoys/data/7mhf/38_7mhf_decoys/merge_2000.csv")
    decoy_df = pd.read_csv(decoy_meta_file, index_col=0)

    pdb_dir = Path("/wynton/group/sali/mhancock/xray/decoys/data/7mhf/38_7mhf_decoys/merge_2000")
    pdb_dir.mkdir()

    for i in range(len(decoy_df)):
        decoy_file = Path(decoy_df["0"].iloc[i])
        new_file = Path(pdb_dir, "{}.pdb".format(i))

        print(decoy_file)
        shutil.copy(
            src=decoy_file,
            dst=new_file
        )

