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
    data_dir = Path(Path.home(), "xray/sample_bench/data/7mhf/179_exp")
    summary_df_file = Path(data_dir, "summary.csv")
    summary_df = pd.read_csv(summary_df_file, index_col=0)

    pdb_dir = Path(Path.home(), "xray/sample_bench/data/7mhf/179_exp/summary")
    f_obs_dir = Path(Path.home(), "xray/data/cifs/7mhf")
    cif_dir = Path(data_dir, "summary_cifs")

    for pdb_file in pdb_dir.glob("*.pdb"):
        cif_name = pdb_file.stem[0:4]
        N = int(pdb_file.stem.split("_")[1][1:])
        J = int(pdb_file.stem.split("_")[2][1:])

        print(pdb_file, cif_name, N, J)

        f_obs_file = Path(f_obs_dir, "{}.cif".format(cif_name))

        pdb_entry = summary_df[(summary_df["cif_name"] == cif_name) & (summary_df["N"] == N) & (summary_df["J"] == J)].iloc[0]

        occs = [pdb_entry["w_{}".format(i)] for i in range(N)]
        print(occs)

        f_model, status_array = get_f_model_from_f_obs(
            pdb_file=pdb_file,
            cif_file=f_obs_file,
            occs=occs
        )

        out_cif_file = Path(cif_dir, "{}.cif".format(pdb_file.stem))
        print(out_cif_file)
        write_cif(
            f_obs=f_model,
            status_array=status_array,
            cif_file=out_cif_file
        )


        # f_obs_file = Path(f_obs_dir, "{}.cif".format(cif_name))
