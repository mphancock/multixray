from pathlib import Path

import IMP
import IMP.atom


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7.pdb")
    s = IMP.atom.AllPDBSelector()
    m = IMP.Model()

    h = IMP.atom.read_pdb(str(pdb_file), s, m)