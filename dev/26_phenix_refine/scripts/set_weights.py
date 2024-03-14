from pathlib import Path
import sys

sys.path.append(str(Path(Path.home(), "xray/src")))
from weights import update_multi_state_model, get_weights_from_pdb_file, get_weights

import IMP
import IMP.atom


if __name__ == "__main__":
    pdb_file = Path("/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/data/n1_n2.pdb")
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    update_multi_state_model(hs=hs, m=m, ws=[0.5502008168734825, 0.4497991831265175])

    # out_pdb_file = Path(str(pdb_file)+".1")
    # print(out_pdb_file)
    IMP.atom.write_multimodel_pdb(hs, str(pdb_file))