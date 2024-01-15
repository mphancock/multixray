from pathlib import Path
import shutil
import math

import IMP
import IMP.atom


if __name__ == "__main__":
    ms = list()
    hs = list()

    m = IMP.Model()

    hs = IMP.atom.read_multimodel_pdb(str(Path(Path.home(), "xray/tmp/2155.pdb")), m, IMP.atom.AllPDBSelector())

    occs = [0.8698016060269730, 0.13019839397302700]

    for i in range(2):
        h = hs[i]

        for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
            at = IMP.atom.Atom(m, pid)
            at.set_occupancy(occs[i])

    IMP.atom.write_multimodel_pdb(hs, str(Path(Path.home(), "xray/tmp/2155_tmp.pdb")))