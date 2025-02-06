import IMP
import IMP.atom

from align_imp import align_one_to_two

"""
Highly specialized function for the refine pipeline for 7MHF where we need to add a ZN atom to the structure for refinement but we don't know what its coordinates and occupancies should be.

We have pdb_file_1 which is a single state and contains a ZN atom. We would like to align the structure from pdb_file_1 to one of the states of pdb_file_2 and then return the new coordinates of the ZN atom and the occupancy.
"""
def get_zn_coords_and_occ_after_align(
    pdb_file_1,
    pdb_file_2
):
    ## add the zn ion if 7mhf series
    m_1, m_2 = IMP.Model(), IMP.Model()
    h_1 = IMP.atom.read_pdb(str(pdb_file_1), m_1, IMP.atom.NonAlternativePDBSelector())
    hs_2 = IMP.atom.read_multimodel_pdb(str(pdb_file_2), m_2, IMP.atom.AllPDBSelector())
    zn_pid_1 = IMP.atom.Selection(h_1, atom_type=IMP.atom.AtomType("HET:ZN  ")).get_selected_particle_indexes()[0]

    ## align the zn atom to the first state of the second model
    align_one_to_two(h_1, hs_2[0])

    zn_coords = IMP.core.XYZ(m_1, zn_pid_1).get_coordinates()
    zn_occ = IMP.atom.Atom(m_1, zn_pid_1).get_occupancy()

    return zn_coords, zn_occ


if __name__ == "__main__":
    from pathlib import Path

    pdb_file_1 = Path(Path.home(), "Documents/xray/data/pdbs/7mhf/7mhf.pdb")
    pdb_file_2 = Path(Path.home(), "Documents/xray/tmp/tmp.pdb")

    zn_coords, zn_occ = get_zn_coords_and_occ_after_align(pdb_file_1=pdb_file_1, pdb_file_2=pdb_file_2)

    print(zn_coords, zn_occ)
