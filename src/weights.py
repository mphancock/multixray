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


def get_weights_from_hs(hs):
    weights = list()
    for h in hs:
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        occ = IMP.atom.Atom(h.get_model(), pids[0]).get_occupancy()
        weights.append(occ)

    return weights


def get_weights(
    floor,
    ws_cur,
    sigma=None

):
    if sigma:
        ws_tmp = np.random.normal(ws_cur, scale=sigma)
    else:
        ws_tmp = [random.random() for _ in range(len(ws_cur))]

    while any(w < floor for w in ws_tmp):
        ws_tmp = [random.random() for _ in range(len(ws_cur))]

    ws_new = [w / sum(ws_tmp) for w in ws_tmp]

    return ws_new