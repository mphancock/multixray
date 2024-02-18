import numpy as np
import random

import IMP
import IMP.core
import IMP.atom


def update_multi_state_model(
        hs,
        m,
        ws
):
    # Update the occupancies of the model based on the weight particle attatched to first pid.
    for i in range(len(hs)):
        occ = ws[i]
        pids = IMP.atom.Selection(hs[i]).get_selected_particle_indexes()
        for pid in pids:
            IMP.atom.Atom(m, pid).set_occupancy(occ)


def get_weights_from_pdb_file(pdb_file):
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    return get_weights_from_hs(hs)


def get_weights_from_hs(hs):
    weights = list()
    for h in hs:
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        occ = IMP.atom.Atom(h.get_model(), pids[0]).get_occupancy()
        weights.append(occ)

    return weights


def get_weights(
    floor,
    n_state,
    occs_cur=None,
    sigma=None

):
    if occs_cur and len(occs_cur) != n_state:
        raise RuntimeError("Length of occs_cur and n_state must be the same")

    if occs_cur:
        occs_tmp = np.random.normal(occs_cur, scale=sigma)
    else:
        occs_tmp = [random.random() for _ in range(n_state)]

    while any(w < floor for w in occs_tmp):
        occs_tmp = [random.random() for _ in range(n_state)]

    ws_new = [w / sum(occs_tmp) for w in occs_tmp]

    return ws_new