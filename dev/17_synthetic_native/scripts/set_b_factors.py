from pathlib import Path

import IMP
import IMP.atom


if __name__ == "__main__":
    for n in [1,2,4]:
        pdb_dir = Path(Path.home(), "xray/dev/19_synthetic_native_2/data/pdbs/{}_state".format(n))

        for pdb_file in pdb_dir.glob("*.pdb"):
            print(pdb_file)
            m = IMP.Model()
            hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.ATOMPDBSelector())

            for h in hs:
                for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
                    at = IMP.atom.Atom(m, pid)
                    at.set_temperature_factor(15)

            IMP.atom.write_multimodel_pdb(hs, str(pdb_file))