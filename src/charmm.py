"""
Author: Matthew Hancock
Date: 08/24/21
Description: Implements a molecular mechanics force field as an IMP restraint that can be used in an IMP sampling/optimization function. The restraint includes bonds, angles, dihedrals, impropers, and non-bonded energy terms using the CHARMM parameterization. The restraint can be used as a restraint class inheriting from IMP.Restraint (CharmmRestraint) or as a list of restraints (get_charmm_restraints).
"""
import sys
from pathlib import Path

import IMP
import IMP.atom
import IMP.core
import IMP.algebra

import align_imp
import miller_ops
import cctbx_score
import update_weights_optimizer_state
from derivatives import get_df_mag_ratio


def get_charmm_restraint_set(m, hs):
    charmm_rs = list()
    for state in range(len(hs)):
        charmm_rs.extend(charmm_restraints(m, hs[state]))

    charmm_r_set = IMP.RestraintSet(m, "CHARMMRestraintSet")
    charmm_r_set.add_restraints(charmm_rs)

    return charmm_r_set


# Return a list of restraints rs on the IMP model m (with hierarchy h) based on potential energy terms from the CHARMM force-field/parameterization.
def charmm_restraints(
        m,
        h
    ):

    rs = list()

    # Configure the IMP model based on the CHARMM parameterization.
    # ff = IMP.atom.get_heavy_atom_CHARMM_parameters()
    # pids = IMP.atom.Selection(h).get_selected_particle_indexes()

    ff = IMP.atom.get_all_atom_CHARMM_parameters()
    topology = ff.create_topology(h)

    ## Why is this not fixing the NTER?
    topology.apply_default_patches()
    # topology.setup_hierarchy(h)
    topology.add_atom_types(h)
    # topology.add_missing_atoms(h)
    # IMP.atom.remove_charmm_untyped_atoms(h)
    # topology.add_coordinates(h)

    ## setup charges first
    charges = topology.add_charges(h)

    # setup waters
    for pid in IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("HET: O  ")).get_selected_particle_indexes():
        charmm_at = IMP.atom.CHARMMAtom.setup_particle(m, pid, "O")
        charge = IMP.atom.Charged.setup_particle(m, pid, -0.834)

    ## for 7mhf
    # ## setup zinc ions
    # for pid in IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("HET:ZN  ")).get_selected_particle_indexes():
    #     IMP.atom.CHARMMAtom.setup_particle(m, pid, "ZN")
    #     charge = IMP.atom.Charged.setup_particle(m, pid, 2.0)

    # add_H_to_N_ter(m, h, topology)

    # for res_id in [41, 80, 163, 172, 246]:
    #     convert_HIS_to_HSP(m, h, ff, topology, res_id)

    bonds = topology.add_bonds(h)
    angles = ff.create_angles(bonds)
    dihedrals = ff.create_dihedrals(bonds)
    impropers = topology.add_impropers(h)

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

    slack = 1
    dist = 3

    # Add a restraint on the non-bonded atoms (Lennard-Jones potential).
    ff.add_radii(h)
    ff.add_well_depths(h)
    atoms = IMP.atom.get_by_type(h, IMP.atom.ATOM_TYPE)

    cont = IMP.container.ListSingletonContainer(m, atoms)
    nbl = IMP.container.ClosePairContainer(cont, dist, slack)

    pair_filter = IMP.atom.StereochemistryPairFilter()
    pair_filter.set_bonds(bonds)
    pair_filter.set_angles(angles)
    pair_filter.set_dihedrals(dihedrals)
    nbl.add_pair_filter(pair_filter)
    sf = IMP.atom.ForceSwitch(6.0, 7.0)
    ljps = IMP.atom.LennardJonesPairScore(sf)
    rs.append(IMP.container.PairsRestraint(ljps, nbl, "nbd"))

    # pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    # r_ljps = PairsRestraintWrapper(pids, ljps, nbl, "nbd")
    # rs.append(r_ljps)

    # sf = IMP.atom.ForceSwitch(6.0, 7.0)
    # cps = IMP.atom.CoulombPairScore(sf)
    # r_cps = IMP.container.PairsRestraint(cps, nbl, "eps")
    # rs.append(r_cps)

    return rs


class CHARMMDerivHolder():
    def __init__(self):
        self.df_dxs = dict()

    def get_df_dxs(self):
        return self.df_dxs

    def get_mean_magnitude(self):
        mean_mag = 0
        for pid in self.df_dxs.keys():
            mean_mag += self.df_dxs[pid].get_magnitude()
        mean_mag /= len(self.df_dxs.keys())

        return mean_mag

    def set_df_dxs(self, df_dxs):
        self.df_dxs = df_dxs




