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
    pdb_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs")
    cif_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/1")
    pdb_files = list(pdb_dir.glob("*.pdb"))
    w_set = 1

    native_df = pd.read_csv(Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/natives.csv"), index_col=0)
    for i in range(len(pdb_files)):
        pdb_file = pdb_files[i]
        model_cif_file = Path(cif_dir, "{}.cif".format(pdb_file.stem))

        w_0 = native_df.loc[int(pdb_file.stem), "weight_{}_0".format(w_set)]
        w_1 = native_df.loc[int(pdb_file.stem), "weight_{}_1".format(w_set)]
        print(w_0, w_1)

        f_obs_array = generate_fmodel.get_f_model(
            pdb_file=pdb_file,
            uc_dim=(58.3050, 36.1540, 25.3620, 90.0000, 103.0900, 90.0000),
            sg_symbol="C 1 2 1",
            res=1,
            ws=(w_0, w_1)
        )

        f_obs_array = generate_fmodel.randomize_amplitude(
            f_obs=f_obs_array
        )

        flags_array = f_obs_array.generate_r_free_flags(
            fraction=.1
        )

        status_array = generate_fmodel.get_status_array(
            flags_array=flags_array
        )

        print(status_array.data())
        print(flags_array.data())

        generate_fmodel.write_cif(
            f_obs=f_obs_array,
            status_array=status_array,
            cif_file=model_cif_file
        )
