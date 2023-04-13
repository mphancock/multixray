"""
Author: Matthew Hancock
Date: 08/24/21
Description: Implements a molecular mechanics force field as an IMP restraint that can be used in an IMP sampling/optimization function. The restraint includes bonds, angles, dihedrals, impropers, and non-bonded energy terms using the CHARMM parameterization. The restraint can be used as a restraint class inheriting from IMP.Restraint (CharmmRestraint) or as a list of restraints (get_charmm_restraints).
"""
import sys
from pathlib import Path
import IMP
import IMP.atom
sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


# Return a list of restraints rs on the IMP model m (with hierarchy h) based on potential energy terms from the CHARMM force-field/parameterization.
def charmm_restraints(
        m,
        h,
        eps=False
):

    rs = list()

    # Configure the IMP model based on the CHARMM parameterization.
    ff = IMP.atom.get_heavy_atom_CHARMM_parameters()
    topology = ff.create_topology(h)

    # topology.apply_default_patches()
    topology.setup_hierarchy(h)

    IMP.atom.remove_charmm_untyped_atoms(h)
    # topology.add_missing_atoms(h)
    topology.add_coordinates(h)
    bonds = topology.add_bonds(h)
    angles = ff.create_angles(bonds)
    dihedrals = ff.create_dihedrals(bonds)
    impropers = topology.add_impropers(h)
    charges = topology.add_charges(h)

    # Add a restraint on the bond lengths.
    cont = IMP.container.ListSingletonContainer(m, bonds, "bnd")
    bss = IMP.atom.BondSingletonScore(IMP.core.Harmonic(0, 1))
    r = IMP.container.SingletonsRestraint(bss, cont, "bnd")
    rs.append(r)

    # Add a restraint on the bond angles.
    cont = IMP.container.ListSingletonContainer(m, angles, "ang")
    bss = IMP.atom.AngleSingletonScore(IMP.core.Harmonic(0, 1))
    r = IMP.container.SingletonsRestraint(bss, cont, "ang")
    rs.append(r)

    # Add a restraint on the dihedral angles.
    cont = IMP.container.ListSingletonContainer(m, dihedrals, "dih")
    bss = IMP.atom.DihedralSingletonScore()
    r = IMP.container.SingletonsRestraint(bss, cont, "dih")
    rs.append(r)

    # Add a restraint on the improper dihedrals (out of plane bending).
    cont = IMP.container.ListSingletonContainer(m, impropers, "imp")
    bss = IMP.atom.ImproperSingletonScore(IMP.core.Harmonic(0, 1))
    rs.append(IMP.container.SingletonsRestraint(bss, cont, "imp"))

    if eps:
        # Add a restraint representing electrostatic forces.
        atoms = IMP.atom.get_by_type(h, IMP.atom.ATOM_TYPE)
        cont = IMP.container.ListSingletonContainer(m, atoms)
        nbl = IMP.container.ClosePairContainer(cont, 10)
        pair_filter = IMP.atom.StereochemistryPairFilter()
        pair_filter.set_bonds(bonds)
        pair_filter.set_angles(angles)
        pair_filter.set_dihedrals(dihedrals)
        nbl.add_pair_filter(pair_filter)
        sf = IMP.atom.ForceSwitch(6.0, 7.0)
        cps = IMP.atom.CoulombPairScore(sf)
        r_cps = IMP.container.PairsRestraint(cps, nbl, "eps")
        rs.append(r_cps)

    # Add a restraint on the non-bonded atoms (Lennard-Jones potential).
    ff.add_radii(h)
    ff.add_well_depths(h)
    atoms = IMP.atom.get_by_type(h, IMP.atom.ATOM_TYPE)
    cont = IMP.container.ListSingletonContainer(m, atoms)
    nbl = IMP.container.ClosePairContainer(cont, 10)
    pair_filter = IMP.atom.StereochemistryPairFilter()
    pair_filter.set_bonds(bonds)
    pair_filter.set_angles(angles)
    pair_filter.set_dihedrals(dihedrals)
    nbl.add_pair_filter(pair_filter)
    sf = IMP.atom.ForceSwitch(6.0, 7.0)
    ljps = IMP.atom.LennardJonesPairScore(sf)
    rs.append(IMP.container.PairsRestraint(ljps, nbl, "nbd"))

    return rs


def ca_distance_restraints(
        m,
        h,
        k
):
    rs = list()

    stiffness = 2

    pids = list()

    cas = IMP.atom.Selection(
        hierarchy=h,
        atom_type=IMP.atom.AT_C
    )
    ns = IMP.atom.Selection(
        hierarchy=h,
        atom_type=IMP.atom.AT_N
    )

    pids.extend(cas.get_selected_particle_indexes())
    pids.extend(ns.get_selected_particle_indexes())

    print(len(cas.get_selected_particle_indexes()))

    for pid in pids:
        xyz = IMP.core.XYZR(m, pid).get_coordinates()
        p = m.get_particle(pid)
        ub = IMP.core.Harmonic(0, k)
        ss = IMP.core.DistanceToSingletonScore(ub, xyz)
        r = IMP.core.SingletonRestraint(m, ss, p)
        rs.append(r)

    return rs


def get_ff_score(
        hs,
        term
):
    score_tot = 0
    m = hs[0].get_model()

    for i in range(len(hs)):
        h = hs[i]
        # w = align_imp.get_pdb_weight(h)
        w = 1 / len(hs)

        ff_scores = dict()
        ff_scores["ff"] = 0

        charmm_rs = charmm_restraints(
            m=m,
            h=h,
            eps=False
        )

        for j in range(len(charmm_rs)):
            term_score = charmm_rs[j].evaluate(
                calc_derivs=False
            )

            ff_scores[charmm_rs[j].get_name()] = term_score
            ff_scores["ff"] = ff_scores["ff"] + term_score

        score = ff_scores[term]

        score_tot = score_tot + score*w

    return score_tot

