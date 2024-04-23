from pathlib import Path
import pandas as pd
import shutil


if __name__ == "__main__":
    cif_df = pd.DataFrame()

    native_df = pd.read_csv(Path(Path.home(), "Documents/xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv"))

    cif_id = 0
    for i in range(10):
        pdb_dir = Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30")
        pdb_file = Path(pdb_dir, "{}.pdb".format(i))

        cif_dir = Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/{}".format(i))
        # cif_dir.mkdir(parents=True, exist_ok=True)

        cif_0 = Path(cif_dir, "0.cif".format(i))
        cif_1 = Path(cif_dir, "1.cif".format(i))

        occs_0_0 = native_df.loc[i, "w_0_0"]
        occs_1_0 = native_df.loc[i, "w_1_0"]
        occs_0_1 = native_df.loc[i, "w_0_1"]
        occs_1_1 = native_df.loc[i, "w_1_1"]

        for j in range(2):
            if j == 0:
                cif_df.loc[cif_id, "cifs"] = str(cif_0)
                cif_df.loc[cif_id, "refs"] = str(pdb_file)
                cif_df.loc[cif_id, "ref_occs"] = str(occs_0_0) + "," + str(occs_1_0)
            elif j == 1:
                cif_df.loc[cif_id, "cifs"] = str(cif_0) + "," + str(cif_1)
                cif_df.loc[cif_id, "refs"] = str(pdb_file) + "," + str(pdb_file)
                cif_df.loc[cif_id, "ref_occs"] = str(occs_0_0) + "," + str(occs_1_0) + ";" + str(occs_0_1) + "," + str(occs_1_1)

            cif_id = cif_id + 1

    cif_df.to_csv(Path(Path.home(), "Documents/xray/dev/29_synthetic_native_3/data/cifs/csvs/7mhf_30.csv"))









