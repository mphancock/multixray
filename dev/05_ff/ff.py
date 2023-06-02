import sys
from pathlib import Path

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/dev/scratch/one_res.pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    for pid in pids:
        IMP.atom.Atom(m, pid).set_occupancy(0)

    rs = charmm_restraints(
        m=m,
        h=h,
        eps=False
    )

    rset_charmm = IMP.RestraintSet(m, 1.0)
    rset_charmm.add_restraints(
        charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
    )

