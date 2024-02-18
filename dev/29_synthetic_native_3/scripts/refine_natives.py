from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(Path.home(), "xray/dev/11_refine/scripts")))
import refine

import IMP
import IMP.atom

if __name__ == "__main__":
    native_df = pd.read_csv(str(Path(Path.home(), "xray/dev/29_synthetic_native_3/data/scores/7mhf_orig.csv")))

    for i in range(len(native_df)):
        pdb_file = native_df.loc[i, "pdb"]
        print(pdb_file)

        out_pdb_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs/7mhf/{}.pdb".format(i))
        n_steps = 10
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        refine.refine(
            hs=hs,
            n_step=n_steps
        )

        IMP.atom.write_multimodel_pdb(hs, str(out_pdb_file))