from pathlib import Path
import numpy as np

import IMP
import IMP.atom

import iotbx.pdb
import mmtbx.model

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


def get_u_anisos_from_file(
    u_aniso_file
):
    # The particles will line up based on the order they are read in from the pdb file.
    model = mmtbx.model.manager(model_input=iotbx.pdb.input(str(u_aniso_file)))
    xray_struct = model.get_xray_structure()
    scatterers = xray_struct.scatterers()

    u_anisos = list()
    for i in range(len(scatterers)):
        scatterer = scatterers[i]
        u_aniso = scatterer.u_star
        u_anisos.append(u_aniso)

    return u_anisos
