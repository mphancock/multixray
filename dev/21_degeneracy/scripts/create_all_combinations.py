from pathlib import Path

import IMP
import IMP.core
import IMP.atom


if __name__ == "__main__":
    data_dir = Path(Path.home(), "xray/dev/21_degeneracy/data")

    pdb_file_1 = Path(data_dir, "AA.pdb")
    pdb_file_2 = Path(data_dir, "output_0/pdbs/90.pdb")

    m_1, m_2 = IMP.Model(), IMP.Model()
    h_1 = IMP.atom.read_pdb(str(pdb_file_1), m_1, IMP.atom.AllPDBSelector())
    h_2 = IMP.atom.read_pdb(str(pdb_file_2), m_2, IMP.atom.AllPDBSelector())

    for first, second in [(1,1),(1,2),(2,1),(2,2)]:
        m = IMP.Model()
        h = IMP.atom.read_pdb(str(pdb_file_1), m, IMP.atom.AllPDBSelector())
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()

        xyzs = list()
        if first == 1:
            pids_ref = IMP.atom.Selection(h_1, residue_index=71).get_selected_particle_indexes()
            xyzs.extend(IMP.core.XYZ(m_1, pid).get_coordinates() for pid in pids_ref)
        else:
            pids_ref = IMP.atom.Selection(h_2, residue_index=71).get_selected_particle_indexes()
            xyzs.extend(IMP.core.XYZ(m_2, pid).get_coordinates() for pid in pids_ref)
        if second == 1:
            pids_ref = IMP.atom.Selection(h_1, residue_index=72).get_selected_particle_indexes()
            xyzs.extend(IMP.core.XYZ(m_1, pid).get_coordinates() for pid in pids_ref)
        else:
            pids_ref = IMP.atom.Selection(h_2, residue_index=72).get_selected_particle_indexes()
            xyzs.extend(IMP.core.XYZ(m_2, pid).get_coordinates() for pid in pids_ref)

        for i in range(len(pids)):
            d = IMP.core.XYZ(m, pids[i]).set_coordinates(xyzs[i])

        for pid in pids:
            IMP.atom.Atom(m, pid).set_occupancy(.5)
            IMP.atom.Atom(m, pid).set_temperature_factor(0)

        out_file = Path(data_dir, "{}{}.pdb".format(first, second))
        IMP.atom.write_pdb(h, str(out_file))

    for pdb_file_names in [["11", "22"],["12", "21"]]:
        hs = list()
        ms = list()
        for pdb_file_name in pdb_file_names:
            pdb_file = Path(data_dir, "{}.pdb".format(pdb_file_name))

            m = IMP.Model()
            h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

            hs.append(h)
            ms.append(m)

        out_file = Path(data_dir, "{}{}.pdb".format(pdb_file_names[0], pdb_file_names[1]))
        IMP.atom.write_multimodel_pdb(hs, str(out_file))






