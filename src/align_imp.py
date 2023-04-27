import IMP
import IMP.atom
import IMP.core
import IMP.algebra
import math
from pathlib import Path
import numpy as np
import sys

import average_structure


# Align pdb file 1 to pdb file 2 and save to pdb file 3.
def align_pdb_files(pdb_1, pdb_2, pdb_3, ca=False):
    sel = IMP.atom.AndPDBSelector(IMP.atom.NonWaterNonHydrogenPDBSelector(), IMP.atom.ATOMPDBSelector())
    # if ca:
    #     # sel = IMP.atom.AndPDBSelector(sel, IMP.atom.ChainPDBSelector("A"))
    #     sel = IMP.atom.AndPDBSelector(sel, IMP.atom.CAlphaPDBSelector())

    m1, m2 = IMP.Model(), IMP.Model()
    h1 = IMP.atom.read_pdb(str(pdb_1), m1, sel)
    h2 = IMP.atom.read_pdb(str(pdb_2), m2, sel)

    try:
        rmsd = align_imp_models(
            h_1=h1,
            m_1=m1,
            ref_h=h2,
            ref_m=m2,
            ca=ca
        )
    except RuntimeError:
        return False

    IMP.atom.write_pdb(h1, str(pdb_3))
    return rmsd


# Align imp molecular hierarchy h to h_ref.
def align(
        h,
        h_0
):
    m = h.get_model()
    m_0 = h_0.get_model()

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    pids_0 = IMP.atom.Selection(h_0).get_selected_particle_indexes()

    xyzs = [IMP.core.XYZR(m, pid) for pid in pids]
    xyzs_0 = [IMP.core.XYZR(m_0, pid) for pid in pids_0]

    if len(xyzs_0) != len(xyzs):
        raise RuntimeError("Length of reference XYZs ({}) and model XYZs ({}) are not equal".format(len(xyzs_ref), len(xyzs)))

    trans_align = IMP.algebra.get_transformation_aligning_first_to_second(
        source=[xyz.get_coordinates() for xyz in xyzs],
        target=[xyz.get_coordinates() for xyz in xyzs_0]
    )

    # Transform all atomic coordinates of h.
    for xyz in xyzs:
        xyz.set_coordinates(trans_align.get_transformed(xyz.get_coordinates()))

    rmsd = IMP.atom.get_rmsd(xyzs, xyzs_0)

    return rmsd


"""
Given 2 sets of input ids of equal length, the function returns a list of all possible matchings between the ids of the first and second input lists. Each matching is a list containing tuples (id_1, id_2).
"""
def get_pairings(
        ids_1,
        ids_2
):
    pairings = list()
    if len(ids_1) == 1:
        pairing = [(ids_1[0], ids_2[0])]
        pairings.append(pairing)
    else:
        for i in range(len(ids_1)):
            id_1 = ids_1[0]
            id_2 = ids_2[i]

            ids_1_tmp = ids_1.copy()
            ids_1_tmp.remove(id_1)
            ids_2_tmp = ids_2.copy()
            ids_2_tmp.remove(id_2)

            tmp_pairings = get_pairings(
                ids_1=ids_1_tmp,
                ids_2=ids_2_tmp
            )

            for tmp_pairing in tmp_pairings:
                pairing = [(id_1,id_2)]
                pairing.extend(tmp_pairing)
                pairings.append(pairing)

    return pairings


"""
Given an IMP hierarchy and model, the function returns the weight of the first atom of each structure in the pdb file.
"""
def get_pdb_weight(
        h
):
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()

    m = h.get_model()
    at = IMP.atom.Atom(m, pids[0])
    occ = at.get_occupancy()

    return occ


"""
Given 2 sets of IMP hierarchies, the function returns the minimum weighted RMSD between the 2 sets of structures. All possible sets of partners between the first and second set are tested. The weighted RMSD is computed as the weighted sum of the matched partner structures based on the weight of the first structure.
"""
def compute_min_weighted_rmsd(
        hs_1,
        hs_2
):
    if len(hs_1) != len(hs_2):
        raise RuntimeError("List of hierarchies not of equal sizes: {} and {}".format(len(hs_1), len(hs_2)))

    ids_1, ids_2 = list(range(len(hs_1))), list(range(len(hs_2)))
    pairings = get_pairings(
        ids_1=ids_1,
        ids_2=ids_2
    )

    min_rmsd = math.inf
    for i in range(len(pairings)):
        total_rmsd = 0
        pairing = pairings[i]
        for j in range(len(pairing)):
            id_1, id_2 = pairing[j]
            h_1 = hs_1[id_1]
            h_2 = hs_2[id_2]
            pair_rmsd = compute_rmsd(
                h=h_1,
                h_0=h_2
            )

            w_1 = get_pdb_weight(h_1)

            total_rmsd = total_rmsd + w_1 * pair_rmsd

        if total_rmsd < min_rmsd:
            min_rmsd = total_rmsd

        # print(total_rmsd)

    return min_rmsd

"""
Given 2 sets of IMP hierarchies, the function returns the RMSD of the average structure of each set of hierarchies. This allows RMSD to be computed between sets of hierarchies of unequal size (eg, a decoy multi-state structure against a single-state native).
"""
def compute_rmsd_between_average(
        hs_1,
        hs_2
):
    avg_dict_1 = average_structure.get_coord_avg_dict(
        hs=hs_1
    )

    avg_dict_2 = average_structure.get_coord_avg_dict(
        hs=hs_2
    )

    rmsd = 0
    for key in avg_dict_1.keys():
        xyz_1 = avg_dict_1[key]
        xyz_2 = avg_dict_2[key]
        mag = (xyz_1-xyz_2).get_magnitude()
        rmsd = rmsd+mag**2

    rmsd = rmsd/len(avg_dict_1.keys())
    rmsd = np.sqrt(rmsd)

    return rmsd

# Compute the root mean squared deviation between 2 IMP models.
def pool_compute_rmsd(
        params
):
    pdb_file_0 = params["pdb_file_0"]
    pdb_file_1 = params["pdb_file_1"]

    m_0, m_1 = IMP.Model(), IMP.Model()
    h_0 = IMP.atom.read_pdb(str(pdb_file_0), m_0, IMP.atom.AllPDBSelector())
    h_1 = IMP.atom.read_pdb(str(pdb_file_1), m_1, IMP.atom.AllPDBSelector())

    rmsd = compute_rmsd(
        h_0=h_0,
        h_1=h_1
    )

    return rmsd, pdb_file_0, pdb_file_1

def compute_rmsd(
        h_0,
        h_1
):
    m_0 = h_0.get_model()
    m_1 = h_1.get_model()

    pids_0 = IMP.atom.Selection(h_0).get_selected_particle_indexes()
    pids_1 = IMP.atom.Selection(h_1).get_selected_particle_indexes()

    xyzs_0 = [IMP.core.XYZR(m_0, pid) for pid in pids_0]
    xyzs_1 = [IMP.core.XYZR(m_1, pid) for pid in pids_1]

    rmsd = IMP.atom.get_rmsd(xyzs_0, xyzs_1)

    return rmsd


def pool_compute_rmsd_aligning_first_to_second(
        params
):
    pdb_file_0 = params["pdb_file_0"]
    pdb_file_1 = params["pdb_file_1"]

    m_0, m_1 = IMP.Model(), IMP.Model()
    h_0 = IMP.atom.read_pdb(str(pdb_file_0), m_0, IMP.atom.AllPDBSelector())
    h_1 = IMP.atom.read_pdb(str(pdb_file_1), m_1, IMP.atom.AllPDBSelector())

    rmsd = compute_rmsd_aligning_first_to_second(
        h_0=h_0,
        h_1=h_1
    )

    return rmsd


def compute_rmsd_aligning_first_to_second(
        h_0,
        h_1
):
    m_0 = h_0.get_model()
    m_1 = h_1.get_model()

    pids_0 = IMP.atom.Selection(h_0).get_selected_particle_indexes()
    pids_1 = IMP.atom.Selection(h_1).get_selected_particle_indexes()

    xyzs_0 = [IMP.core.XYZR(m_0, pid) for pid in pids_0]
    xyzs_1 = [IMP.core.XYZR(m_1, pid) for pid in pids_1]

    if len(xyzs_0) != len(xyzs_1):
        raise RuntimeError("Length of reference XYZs ({}) and model XYZs ({}) are not equal".format(len(xyzs_0), len(xyzs_1)))

    t_align = IMP.algebra.get_transformation_aligning_first_to_second(
        source=[xyz.get_coordinates() for xyz in xyzs_0],
        target=[xyz.get_coordinates() for xyz in xyzs_1]
    )
    rmsd = IMP.atom.get_rmsd_transforming_first(t_align, xyzs_0, xyzs_1)

    return rmsd


if __name__ == "__main__":
    pdb_file_ref = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    pdb_file = Path(Path.home(), "xray/tmp/66_0_1000.pdb")

    m_0, m = IMP.Model(), IMP.Model()
    h_0 = IMP.atom.read_pdb(str(pdb_file_ref), m_0, IMP.atom.AllPDBSelector())
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    rmsd = align(
        h=h,
        h_0=h_0
    )
    print(rmsd)

    IMP.atom.write_pdb(h, str(Path(Path.home(), "xray/tmp/66_0_1000_align.pdb")))

    #
    # rmsd = compute_rmsd(
    #     h_1=h_0,
    #     h_2=h_1
    # )
    # print(rmsd)
