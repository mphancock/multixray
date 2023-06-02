from pathlib import Path
import sys

import IMP
import IMP.atom


if __name__ == "__main__":
    ## READ IN WATER.
    pdb_file = Path(Path.home(), "xray/dev/01_h20/output_0/pdbs/0.pdb")
    s = IMP.atom.AllPDBSelector()
    m = IMP.Model()

    h = IMP.atom.read_multimodel_pdb(str(pdb_file), m, s)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    print(len(pids))
