from pathlib import Path
import sys
import numpy as np
import pandas as pd

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import cctbx.miller
import iotbx
from scitbx.array_family import flex

sys.path.append(str(Path(Path.home(), "xray/data/cifs/scripts")))
import generate_fmodel


if __name__ == "__main__":
    n_state = 2
    n_cond = 2

    native_df = pd.read_csv(Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/7mhf_20.csv"), index_col=0)

    status_array=None

    for i in range(len(native_df)):
        cif_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/7mhf_20_raw/{}".format(i))
        cif_dir.mkdir(exist_ok=True)

        for cond in list(range(n_cond)):
            pdb_file = native_df.loc[i, "pdb"]
            print(pdb_file)

            model_cif_file = Path(cif_dir, "{}.cif".format(cond))

            occs = list()
            for state in list(range(n_state)):
                occs.append(native_df.loc[i, "w_{}_{}".format(state, cond)])

            print(occs)

            f_obs_array = generate_fmodel.get_f_model(
                pdb_file=pdb_file,
                uc_dim=(114.968, 54.622, 45.194, 90.000, 101.675, 90.000),
                sg_symbol="C 1 2 1",
                res=1,
                ws=occs
            )

            # f_obs_array = generate_fmodel.randomize_amplitude(
            #     f_obs=f_obs_array,
            #     dist="uni",
            #     std=100
            # )

            if not status_array:
                flags_array = f_obs_array.generate_r_free_flags(
                    fraction=.05,
                    max_free=20000
                )

                status_array = generate_fmodel.get_status_array(
                    flags_array=flags_array
                )

            # print(status_array.data())
            # print(flags_array.data())

            generate_fmodel.write_cif(
                f_obs=f_obs_array,
                status_array=status_array,
                cif_file=model_cif_file
            )
