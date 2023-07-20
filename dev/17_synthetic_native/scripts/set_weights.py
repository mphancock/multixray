from pathlib import Path

import IMP
import IMP.atom


if __name__ == "__main__":
    pdb_dir = Path(Path.home(), "xray/dev/17_synthetic_native/data/pdbs/2_state_ref")
    new_pdb_dir = Path(Path.home(), "xray/dev/17_synthetic_native/data/pdbs/2_state_uni")

    for pdb_file in pdb_dir.glob("*.pdb"):
        print(pdb_file)
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        for h in hs:
            pids = IMP.atom.Selection(h).get_selected_particle_indexes()
            for pid in pids:
                IMP.atom.Atom(m, pid).set_occupancy(1/len(hs))

        IMP.atom.write_multimodel_pdb(hs, str(Path(new_pdb_dir, pdb_file.name)))