import IMP
import IMP.atom


def position_restraints(
        # h_ref,
        # m_ref,
        m,
        h,
        k
):
    # pids_ref = IMP.atom.Selection(h_ref).get_selected_particle_indexes()
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    rs = list()

    for pid in pids:
        xyz = IMP.core.XYZR(m, pid).get_coordinates()
        p = m.get_particle(pid)
        ub = IMP.core.Harmonic(0, k)
        ss = IMP.core.DistanceToSingletonScore(ub, xyz)
        r = IMP.core.SingletonRestraint(m, ss, p)
        rs.append(r)

    return rs