from pathlib import Path
import sys

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "Documents/xray/src")))
import align_imp


if __name__ == "__main__":
    pdb_dir = Path(Path.home(), "Documents/xray/dev/17_synthetic_native/data/pdbs/2_state_ref")

    for pdb_file in  pdb_dir.glob("*"):
        order_pdb_file = Path(Path.home(), "Documents/xray/dev/17_synthetic_native/data/pdbs/2_state_ref_ordered", pdb_file.name)
        print(pdb_file)

        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        ordered_hs = align_imp.get_ordered_hs(hs)

        IMP.atom.write_multimodel_pdb(ordered_hs, str(order_pdb_file))




