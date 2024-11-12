from pathlib import Path

import numpy as np

import IMP
import IMP.atom

import iotbx.pdb


def get_string_from_sample_sched(sa_sched):
    sa_str = ""
    for sa_step in sa_sched:
        step_str = ""
        for key, val in sa_step.items():
            if key == "res":
                val_str = "{:.2f}".format(val)  # Format float values with 2 decimal points
            else:
                val_str = str(val)
            step_str += f"{key}{val_str},"
        sa_str += step_str[:-1] + ";"
    return sa_str[:-1]


def get_n_state_from_pdb_file(pdb_file):
    pdb_return = pool_read_pdb(pdb_file)

    if isinstance(pdb_return, RuntimeError):
        raise pdb_return
    else:
        m, hs = pdb_return

    return len(hs)


def get_residue_indexes(h): # h is a Hierarchy
    res_ids = list()
    for pid in res_ids:
        if IMP.atom.Residue.get_is_setup(h.get_particle(pid)):
            res = IMP.atom.Residue(h.get_particle(pid))
            res_ids.append(res.get_index())

    res_ids.sort()

    return res_ids


def pool_read_pdb(
    pdb_file
):
    # print(pdb_file)

    if pdb_file.exists():
        m = IMP.Model()

        try:
            hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        except IMP.ValueException:
            return RuntimeError("{} cannot be read".format(pdb_file))

        return m, hs
    else:
        return RuntimeError("{} file does not exist".format(pdb_file))


def get_crystal_info(pdb_file):
    pdb_input = iotbx.pdb.hierarchy.input(file_name=str(pdb_file))
    structure = pdb_input.xray_structure_simple()
    symmetry = structure.crystal_symmetry()
    uc = symmetry.unit_cell()
    uc_dimensions = uc.parameters()
    sg_info = symmetry.space_group_info()
    sg_symbol = sg_info.symbol_and_number().split("(")[0]

    return uc_dimensions, sg_symbol


def get_ordered_hs(
        h_0s,
        occs
):
    # Need to invert the list.
    ids = list(np.argsort(occs))
    ids.reverse()
    h_0s_ordered = [h_0s[i] for i in ids]
    return h_0s_ordered
