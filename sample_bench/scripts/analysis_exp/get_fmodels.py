from pathlib import Path
import pandas as pd
import numpy as np
import sys
import multiprocessing

import mmtbx.f_model
import mmtbx.model
import cctbx.crystal
import cctbx.xray
import iotbx

sys.path.append(str(Path(Path.home(), "xray/src")))
# sys.path.append("../src")
import miller_ops
sys.path.append(str(Path(Path.home(), "xray/data/cifs/scripts")))
from generate_fmodel import get_f_model_from_f_obs, write_cif


if __name__ == "__main__":
    data_dir = Path(Path.home(), "xray/sample_bench/data/7mhf/189_exp_ref_10000")
    summary_df_file = Path(data_dir, "summary_best.csv")
    summary_df = pd.read_csv(summary_df_file, index_col=0)

    # pdb_dir = Path(Path.home(), "xray/sample_bench/data/7mhf/179_exp/summary")
    pdb_df = pd.read_csv(Path(data_dir, "summary_best.csv"), index_col=0)
    f_obs_dir = Path(Path.home(), "xray/data/cifs/7mhf")
    f_model_dir = Path(data_dir, "summary_best_cif")

    for index in pdb_df.index:
        N = pdb_df.loc[index, "N"]
        J = pdb_df.loc[index, "J"]

        cif_name = None
        for cif_name in ["7mhf", "7mhg", "7mhh", "7mhi", "7mhj", "7mhk"]:
            if not np.isnan(pdb_df.loc[index, "xray_{}".format(cif_name)]):
                break


        f_obs_file = Path(f_obs_dir, "{}.cif".format(cif_name))
        f_model_file = Path(f_model_dir, "{}.cif".format(index))

        pdb_file = pdb_df.loc[index, "pdb"]
        # pdb_entry = summary_df[(summary_df["cif_name"] == cif_name) & (summary_df["N"] == N) & (summary_df["J"] == J)].iloc[0]

        occs = [pdb_df.loc[index, "w_{}_{}".format(i, cif_name)] for i in range(N)]
        print(pdb_file, cif_name, N, J)
        print(occs)

        f_model, status_array = get_f_model_from_f_obs(
            pdb_file=pdb_file,
            cif_file=f_obs_file,
            occs=occs
        )

        # out_cif_file = Path(cif_dir, "{}.cif".format(index.stem))
        # print(out_cif_file)
        write_cif(
            f_obs=f_model,
            status_array=status_array,
            cif_file=f_model_file
        )


        # f_obs_file = Path(f_obs_dir, "{}.cif".format(cif_name))
