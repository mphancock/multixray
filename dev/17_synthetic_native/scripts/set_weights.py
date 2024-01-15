from pathlib import Path
import random
import numpy as np

import IMP
import IMP.atom

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
import weights


if __name__ == "__main__":
    pdb_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs/2_state_0")
    new_pdb_dir = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs")

    for pdb_file in pdb_dir.glob("*.pdb"):
        print(pdb_file)
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        ws_cur = weights.get_weights_from_hs(hs)
        ws_new = weights.get_weights(
            floor=.05,
            ws_cur=ws_cur,
            sigma=.05
        )

        ws_new = [1.00, 1.00]
        print(ws_cur, ws_new)

        weights.update_multi_state_model(
            hs=hs,
            m=m,
            ws=ws_new
        )

        IMP.atom.write_multimodel_pdb(hs, str(Path(new_pdb_dir, pdb_file.name)))