from pathlib import Path
import numpy as np

import IMP
import IMP.atom


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