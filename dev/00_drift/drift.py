from pathlib import Path

import IMP
import IMP.atom
import IMP.core


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    drift_pdb_file = Path(Path.home(), "xray/dev/drift/3ca7_clean_drift_1.pdb")

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    for pid in pids:
        xyz = IMP.core.XYZ(m, pid)
        x_0 = xyz.get_coordinate(0)
        xyz.set_coordinate(0, x_0+10)

    IMP.atom.write_pdb(h, str(drift_pdb_file))


