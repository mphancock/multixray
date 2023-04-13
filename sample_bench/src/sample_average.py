from pathlib import Path
import math
import sys

import IMP
import IMP.atom
import IMP.core
import IMP.algebra

sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


def get_coord_avg_dict(
        hs
):
    sel = IMP.atom.Selection(hs[0])
    pids = sel.get_selected_particle_indexes()

    norm_fact = 0
    # We are making a big assumption that the occupancy of all atoms in the structure is the same.
    for h in hs:
        occ = IMP.atom.Atom(h.get_model(), IMP.atom.Selection(h).get_selected_particle_indexes()[0]).get_occupancy()

        norm_fact = norm_fact+occ

    avg_coord_dict = dict()
    for pid in pids:
        coord_avg = IMP.algebra.Vector3D(0,0,0)
        for i in range(len(hs)):
            h = hs[i]
            m = h.get_model()
            at = IMP.atom.Atom(m, pid)

            d = IMP.core.XYZ(m, pid)

            coord_avg = coord_avg + d.get_coordinates()*at.get_occupancy()

        coord_avg = coord_avg / norm_fact
        avg_coord_dict[pid] = coord_avg

    return avg_coord_dict


def get_average_pdb_file_from_pdb_files(
        pdb_files
):
    hs = list()
    ms = list()

    for pdb_file in pdb_files:
        s = IMP.atom.AllPDBSelector()
        m = IMP.Model()
        h = IMP.atom.read_pdb(str(pdb_file), m, s)

        hs.append(h)
        ms.append(m)

    # First, compute the average coordinates.
    coord_avg_dict = get_coord_avg_dict(
        hs=hs
    )

    # Create an artificial structure with the computed average coordinates.
    m_true_avg = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h_synth_avg = IMP.atom.read_pdb(str(pdb_files[0]), m_true_avg, s)

    for pid in coord_avg_dict.keys():
        avg_coords = coord_avg_dict[pid]

        xyz = IMP.core.XYZ(m_true_avg, pid)
        xyz.set_coordinates(avg_coords)

    # Find the structure that has minimal rmsd with the synthetic average structure.
    min_rmsd = math.inf
    min_id = -1
    for i in range(len(hs)):
        rmsd = align_imp.compute_rmsd(
            h=hs[i],
            h_0=h_synth_avg,
            align=False
        )
        print(pdb_files[i], rmsd)

        if rmsd < min_rmsd:
            min_id = i
            min_rmsd = rmsd

    return pdb_files[min_id]


# def get_run_average_pdb_file_dict(
#         pdb_files
# ):
#     run_pdb_file_dict = dict()
#     for pdb_file in pdb_files:
#         run = pdb_file.




if __name__ == "__main__":
    pdb_files = list()
    for i in range(100):
        pdb_files.append(Path(Path.home(), "xray/decoys/data/decoy_sets/3ca7/3ca7_N_1000_x1/{}.pdb".format(i)))
    pdb_files.append(Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb"))

    avg_pdb_file = get_average_pdb_file_from_pdb_files(
        pdb_files=pdb_files
    )
    print(avg_pdb_file)