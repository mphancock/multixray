import sys
from pathlib import Path

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhk/7mhk_clean_h20.pdb")

    m = IMP.Model()
    sel = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, sel)

    charmm_rs = charmm.charmm_restraints(m, h)

    for charmm_r in charmm_rs:
        print(charmm_r.get_name(), charmm_r.evaluate(False))

    for i in range(677):
        res_pids = IMP.atom.Selection(h, residue_index=i).get_selected_particle_indexes()

        atom_types = list()
        for pids in res_pids:
            atom_types.append(IMP.atom.Atom(m, pids).get_atom_type())

        for atom_type in atom_types:
            atom_types_copy = atom_types.copy()
            atom_types_copy.remove(atom_type)
            if atom_type in atom_types_copy:
                print(i, atom_type)


