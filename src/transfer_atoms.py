from pathlib import Path

import IMP
import IMP.atom


"""
copy atom indexed by pid from h_1 to h_2.
"""
def copy_atom(
    pid_1,
    h_1,
    h_2,
):
    m_1 = h_1.get_model()
    m_2 = h_2.get_model()

    at_1 = IMP.atom.Atom(m_1, pid_1)
    res_1 = at_1.get_parent()

    chain_2 = h_2.get_children()[0]
    p_at_1 = IMP.Particle(m_2)
    at_2 = IMP.atom.Atom.setup_particle(p_at_1, at_1)
    at_2.set_temperature_factor(at_1.get_temperature_factor())
    at_2.set_occupancy(at_1.get_occupancy())

    pid_2 = p_at_1.get_index()
    xyz_1 = IMP.core.XYZ(m_1, pid_1)
    xyz_2 = IMP.core.XYZ.setup_particle(m_2, pid_2, xyz_1.get_coordinates())

    res_1 = IMP.atom.get_residue(at_1)
    p_res_2 = IMP.Particle(m_2)
    res_2 = IMP.atom.Residue.setup_particle(p_res_2, res_1)

    res_2.add_child(at_2)
    chain_2.add_child(res_2)


if __name__ == "__main__":
    m_1, m_2 = IMP.Model(), IMP.Model()
    hs_1 = IMP.atom.read_multimodel_pdb(str(Path(Path.home(), "xray/tmp/7mhh_tmp.pdb")), m_1, IMP.atom.AllPDBSelector())
    hs_2 = IMP.atom.read_multimodel_pdb(str(Path(Path.home(), "xray/tmp/500.pdb")), m_2, IMP.atom.AllPDBSelector())

    pids = IMP.atom.Selection(hs_1[0], atom_type=IMP.atom.AtomType("HET:ZN  ")).get_selected_particle_indexes()
    pid = pids[0]

    for h_2 in hs_2:
        align_and_add_atom(pid, hs_1[0], h_2)

    IMP.atom.write_multimodel_pdb(hs_2, str(Path(Path.home(), "xray/tmp/tmp_zn.pdb")))
