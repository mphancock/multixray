from pathlib import Path
import sys

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import average_structure


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/dev/21_degeneracy/data/1122.pdb")
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    coord_avg_dict = average_structure.get_coord_avg_dict(
        hs=hs
    )

    # print(coord_avg_dict)

    pdb_files = Path(Path.home(), "xray/dev/21_degeneracy/data/output_0/pdbs").glob("*.pdb")

    for pdb_file in pdb_files:
        if pdb_file.stem in ["0", "90"]:
            continue

        m_1 = IMP.Model()
        h_1 = IMP.atom.read_pdb(str(pdb_file), m_1, IMP.atom.AllPDBSelector())

        m_2 = IMP.Model()
        h_2 = IMP.atom.read_pdb(str(pdb_file), m_2, IMP.atom.AllPDBSelector())

        pids = IMP.atom.Selection(h_1).get_selected_particle_indexes()
        for pid in pids:
            xyz = IMP.core.XYZ(m_1, pid).get_coordinates()
            xyz_avg = coord_avg_dict[pid]

            xyz_new = list()
            for i in range(3):
                delta = xyz[i] - xyz_avg[i]
                xyz_new.append(xyz_avg[i] - delta)

            # print(xyz, xyz_avg, xyz_new)

            IMP.core.XYZ(m_2, pid).set_coordinates(xyz_new)

            IMP.atom.Atom(m_1, pid).set_occupancy(.5)
            IMP.atom.Atom(m_1, pid).set_temperature_factor(0)
            IMP.atom.Atom(m_2, pid).set_occupancy(.5)
            IMP.atom.Atom(m_2, pid).set_temperature_factor(0)

        degen_pdb_file = Path(Path.home(), "xray/dev/21_degeneracy/data/degens/{}.pdb".format(pdb_file.stem))
        IMP.atom.write_multimodel_pdb([h_1, h_2], str(degen_pdb_file))


