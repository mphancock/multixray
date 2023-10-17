from pathlib import Path
import random

import IMP
import IMP.atom


if __name__ == "__main__":
    pdb_dir = Path(Path.home(), "xray/dev/19_synthetic_native_2/data/pdbs/4_state")
    new_pdb_dir = Path(Path.home(), "xray/dev/19_synthetic_native_2/data/pdbs/4_state_2")

    for pdb_file in pdb_dir.glob("*.pdb"):
        print(pdb_file)
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        weights = list()
        for i in range(len(hs)):
            weights.append(random.random())

        weights = [w/sum(weights) for w in weights]

        for i in range(len(hs)):
            h = hs[i]

            pids = IMP.atom.Selection(h).get_selected_particle_indexes()
            for pid in pids:
                IMP.atom.Atom(m, pid).set_occupancy(weights[i])

        IMP.atom.write_multimodel_pdb(hs, str(Path(new_pdb_dir, pdb_file.name)))