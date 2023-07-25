from pathlib import Path
import random

import IMP
import IMP.atom


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_refine.pdb")

    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()

    h_clone = IMP.atom.create_clone(h)

    pids_clone = IMP.atom.Selection(h_clone)
    IMP.atom.Atom(h_clone.get_model(), pids_clone.get_selected_particle_indexes()[0]).set_occupancy(0.5)

    print(IMP.atom.Atom(h.get_model(), pids[0]).get_occupancy())

    print(h.get_model())
    print(h_clone.get_model())

    if h.get_model() == h_clone.get_model():
        print("True")


