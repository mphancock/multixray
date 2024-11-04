from pathlib import Path
import math
import sys
import pandas as pd
import multiprocessing

import IMP
import IMP.atom
import IMP.core
import IMP.algebra

sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


def get_coords_and_occs(
        pdb_file
):
    s = IMP.atom.AllPDBSelector()
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    sel = IMP.atom.Selection(h)
    pids = sel.get_selected_particle_indexes()

    coords_dict = dict()
    occs_dict = dict()
    for pid in pids:
        d = IMP.core.XYZ(m, pid)
        at = IMP.atom.Atom(m, pid)
        coords_dict[pid] = d.get_coordinates()
        occs_dict[pid] = at.get_occupancy()

    return coords_dict, occs_dict

"""
It is useful to have a version of get_coord_avg_dict where the reading of pdb files (generally the rate limiting step) is multiprocessed. This was implemented for the trajectory analysis case were we are attempting to compute averages from a large number of structures.
"""
def get_coord_avg_dict_from_pdb_files(
        pdb_files
):
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_files[0]), m, IMP.atom.AllPDBSelector())
    sel = IMP.atom.Selection(h)
    pids = sel.get_selected_particle_indexes()

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(
        get_coords_and_occs,
        pdb_files
    )

    coords_dicts, occs_dicts = list(), list()
    for coords_dict, occs_dict in pool_results:
        coords_dicts.append(coords_dict)
        occs_dicts.append(occs_dict)

    # Get normalization factor from the first atom of each structure.
    norm_fact = 0
    for occ_dict in occs_dicts:
        norm_fact = norm_fact + occ_dict[pids[0]]

    # Compute the average coordinate.
    avg_coord_dict = dict()
    for pid in pids:
        coord_avg = IMP.algebra.Vector3D(0,0,0)
        for i in range(len(coords_dicts)):
            coords_dict = coords_dicts[i]
            occs_dict = occs_dicts[i]

            coord_avg = coord_avg + coords_dict[pid]*occs_dict[pid]

        coord_avg = coord_avg / norm_fact
        avg_coord_dict[pid] = coord_avg

    return avg_coord_dict


"""
This function is useful for computing the average structure from a set of structures (hs). It returns a dictionary containing the average coordinates for each atom in the structure.

We check that each h in hs has equal size. Additionally, we do not assume that the occupancies of the structures are normalized or that all occupancies within a structure are equal.

**********
Params
    hs: list of hierarchies to compute the average structure from.

**********
Returns
    avg_coord_dict: dictionary containing the average coordinates for each atom in the structure. The dictionary is indexed by the particle id of the first structure in the set and the value is the average coordinate as a 3D IMP algebra vector.

"""
def get_coord_avg_dict(
        hs,
        occs
):
    if occs and len(occs) != len(hs):
        raise RuntimeError("The length of the weights list is not equal to the number of structures: {} and {}".format(len(occs), len(hs)))

    avg_dict = dict()
    norm_dict = dict()

    pids_0 = IMP.atom.Selection(hs[0]).get_selected_particle_indexes()
    for pid_0 in pids_0:
        avg_dict[pid_0] = IMP.algebra.Vector3D(0,0,0)
        norm_dict[pid_0] = 0

    # Check that all structures are of equal size.
    for h in hs:
        n_pid = len(IMP.atom.Selection(h).get_selected_particle_indexes())
        if n_pid != len(pids_0):
            raise RuntimeError("Structures are not of equal size {} and {}.".format(n_pid, len(pids_0)))

    # for h in hs:
    for i in range(len(hs)):
        h = hs[i]
        m = h.get_model()

        pids = IMP.atom.Selection(h).get_selected_particle_indexes()

        for j in range(len(pids_0)):
            # pid_0 is the reference pid from the first structure that is used to index the dictionary.
            pid_0 = pids_0[j]
            pid = pids[j]

            at = IMP.atom.Atom(m, pid)
            d = IMP.core.XYZ(m, pid)

            # If a weight is provided use it, otherwise use the occupancy of the atom.
            avg_dict[pid_0] = avg_dict[pid_0] + d.get_coordinates()*occs[i]
            norm_dict[pid_0] = norm_dict[pid_0] + occs[i]

    # Normalize the averages.
    for pid_0 in pids_0:
        avg_dict[pid_0] = avg_dict[pid_0] / norm_dict[pid_0]

    return avg_dict


def get_average_pdb_file_from_pdb_files(
        pdb_files
):
    hs = list()
    ms = list()

    coord_avg_dict = get_coord_avg_dict_from_pdb_files(
        pdb_files=pdb_files
    )

    # Create an artificial structure with the computed average coordinates.
    m_synth_avg = IMP.Model()
    h_synth_avg = IMP.atom.read_pdb(str(pdb_files[0]), m_synth_avg, IMP.atom.AllPDBSelector())
    synth_avg_pdb_file = str(Path(Path.home(), "xray/tmp/synth_avg.pdb"))

    for pid in coord_avg_dict.keys():
        avg_coords = coord_avg_dict[pid]

        xyz = IMP.core.XYZ(m_synth_avg, pid)
        xyz.set_coordinates(avg_coords)
    IMP.atom.write_pdb(h_synth_avg, synth_avg_pdb_file)

    # print("BREAK")

    # Find the structure that has minimal rmsd with the synthetic average structure.
    pool_params = list()
    for pdb_file in pdb_files:
        params_dict = dict()
        params_dict["pdb_file_0"] = pdb_file
        params_dict["pdb_file_1"] = synth_avg_pdb_file
        pool_params.append(params_dict)

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(
        align_imp.pool_compute_rmsd,
        pool_params
    )

    rmsd_df = pd.DataFrame(index=pdb_files, columns=["rmsd"], dtype=float)
    for rmsd, pdb_file_0, pdb_file_1 in pool_results:
        rmsd_df.loc[pdb_file_0, "rmsd"] = rmsd

    # print(rmsd_df.head())
    return rmsd_df.idxmin()["rmsd"]


# if __name__ == "__main__":
#     data_dir = Path(Path.home(), "xray/sample_bench/unit_tests/data")
#     pdb_files=[Path(data_dir, "0.pdb"), Path(data_dir, "1.pdb"), Path(data_dir, "2.pdb")]
#     coord_avg_dict = get_coord_avg_dict_from_pdb_files(
#         pdb_files=pdb_files
#     )

#     print(coord_avg_dict.values())

#     avg_pdb_file = get_average_pdb_file_from_pdb_files(
#         pdb_files=pdb_files,
#     )
#     print(avg_pdb_file)

#     # print(avg_rmsd_df.head())