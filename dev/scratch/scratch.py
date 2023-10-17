from pathlib import Path
import shutil
import math

import IMP
import IMP.atom


if __name__ == "__main__":
    ms = list()
    hs = list()

    m = IMP.Model()

    for i in range(2):
        h = IMP.atom.read_pdb(str(Path(Path.home(), "xray/data/pdbs/7mhf/7mhf_refine.pdb")), m)

        hs.append(h)
        ms.append(m)

    pids_1 = IMP.atom.Selection(hs[0]).get_selected_particle_indexes()
    pids_2 = IMP.atom.Selection(hs[1]).get_selected_particle_indexes()

    print(pids_1[0] == pids_2[0])
