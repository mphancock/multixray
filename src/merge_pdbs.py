import IMP
import IMP.atom
from pathlib import Path
import numpy as np
import multiprocessing
import pandas as pd


def write_merge_pdf_file_pool(
        param_dict
):
    write_merge_pdb_file(
        merge_pdb_file=param_dict["merge_pdb_file"],
        pdb_files=param_dict["pdb_files"],
        occs=param_dict["occs"],
        n=param_dict["n"],
        order=param_dict["order"],
        state=param_dict["state"]
    )

    return param_dict["id"]


def read_pdb_pool(
        pdb_file
):
    print(pdb_file)
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    return m, hs


"""
When writing multi-state MD trajectories, we need each model to represent a distinct frame. Thus, we need to relabel each structure in the multi-state model as a distinct chain in the same model.
"""
def relabel_chains(
        hs,
        m
):
    h_0 = hs[0]
    for i in range(1, len(hs)):
        h_0.add_child(IMP.atom.get_root(hs[i]).get_children()[0])

    hchains = list()
    for pid in m.get_particle_indexes():
        if IMP.atom.Chain.get_is_setup(m, pid):
            hchain = IMP.atom.Chain(m, pid)
            # print(hchain.get_id())
            hchains.append(hchain)

    for j in range(len(hchains)):
        hchain = hchains[j]
        chain_id = chr(ord('A') + j)
        hchain.set_id(chain_id)
        # print(hchain.get_id())

    return hs


"""
This function is used to merge multiple pdb files together in a multi-model pdb file. This is useful in creating MD trajectory movies from sampling or merging together multiple 1-state structures into a multi-state structure.

Params:
    merge_pdb_file: Path to the output multi-model pdb file.

    pdb_files: List of paths to the pdb files to merge together.

    occs: List of occupancy values for each pdb file. The length of this list must be equal to the length of pdb_files.

    n: Number of pdb files to merge together. If n < 0, then all pdb files will be merged together. If n > len(pdb_files), then all pdb files will be merged together.

    order: If True, then the pdb files will be merged together in order of their file names. If False, then the pdb files will be merged together in the order they are given in the pdb_files list.

    state: The state of the model to read from the pdb files.

Returns:
    None

"""
def write_merge_pdb_file(
        merge_pdb_file,
        pdb_files,
        occs,
        n,
        order,
        state
):
    hs_all = list()
    ms_all = list()


    # Sort based on the numeric part of the filename
    sorted_pdb_files = sorted(pdb_files, key=lambda x: int(x.stem))

    print(sorted_pdb_files)
    # if order:
    #     pdb_ids = [int(pdb_file.stem) for pdb_file in pdb_files]
    #     pdb_min = np.min(pdb_ids)
    #     pdb_max = np.max(pdb_ids)

    #     if n < 0 or n > pdb_max:
    #         max_it = pdb_max
    #     else:
    #         max_it = n

    #     pdb_dir = pdb_files[0].parents[0]
    #     for i in range(pdb_min, max_it):
    #         pdb_files_order.append(Path(pdb_dir, "{}.pdb".format(i)))
    # else:
    #     pdb_files_order = pdb_files

    # pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    # pool_results = pool_obj.imap(read_pdb_pool, pdb_files_order)

    # i = 0
    # for pool_result in pool_results:
    #     m, hs = pool_result

    #     if state < 0:
    #         hs = relabel_chains(hs, m)
    #         h = hs[0]
    #     else:
    #         h = hs[state]

    #     pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    #     if occs:
    #         for pid in pids:
    #             at = IMP.atom.Atom(m, pid)
    #             at.set_occupancy(occs[i])

    #     hs_all.append(h)
    #     ms_all.append(m)
    #     i = i+1

    hs, ms = list(), list()
    for pdb_file in sorted_pdb_files:
        print(pdb_file)
        m = IMP.Model()
        h = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.NonWaterPDBSelector())[0]
        # h = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())[0]
        ms.append(m)
        hs.append(h)

    # Cannot return hs, because the hierarchy objects will get deallocated.
    IMP.atom.write_multimodel_pdb(hs, str(merge_pdb_file))
