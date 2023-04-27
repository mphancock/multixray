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


if __name__ == "__main__":
    data_dir = Path(Path.home(), "xray/sample_bench/unit_tests/data")
    pdb_files=[Path(data_dir, "0.pdb"), Path(data_dir, "1.pdb"), Path(data_dir, "2.pdb")]
    coord_avg_dict = get_coord_avg_dict_from_pdb_files(
        pdb_files=pdb_files
    )

    print(coord_avg_dict.values())

    avg_pdb_file = get_average_pdb_file_from_pdb_files(
        pdb_files=pdb_files,
    )
    print(avg_pdb_file)

    # print(avg_rmsd_df.head())