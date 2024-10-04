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


import IMP
import IMP.core
import IMP.algebra

import miller_ops
import cctbx_score
import update_weights_optimizer_state
from derivatives import get_df_mag_ratio


## we can't create a wrapper around restraints (cant call evaluate or do_score_add_derviatives) so need to create a fake restraint to track the derivatives
class CHARMMRestraint(IMP.Restraint):
    def __init__(
            self,
            msmc_m
    ):
        IMP.Restraint.__init__(self, msmc_m.get_m(), "CHARMMShadowRestraint%1%")
        self.msmc_m = msmc_m
        self.hs = msmc_m.get_hs()
        self.n_state = len(self.hs)
        self.pids = msmc_m.get_pids()

        # Gradients and scores
        self.df_dxs = dict()
        for pid in self.pids:
            self.df_dxs[pid] = IMP.algebra.Vector3D(0,0,0)
        self.score = 0

        # setup the charmm restraints
        self.charmm_rs = list()
        for state in range(self.n_state):
            self.charmm_rs.extend(charmm_restraints(
                m=self.get_model(),
                h=self.hs[state],
                eps=False
            ))


    def get_df_dict(self):
        return self.df_dxs

    def get_f(self):
        return self.score

    def do_add_score_and_derivatives(self, sa):
        self.score = 0
        for r in self.charmm_rs:
            r.add_score_and_derivatives(sa)
            self.score += r.get_last_score()

        ## create dictionary for the derivatives
        ## derivatives are stored in the order of atoms/scatterers
        # print(sa.get_score())
        for pid in self.pids:
            d = IMP.core.XYZR(self.get_model(), pid)
            self.df_dxs[pid] = d.get_derivatives()

    def do_get_inputs(self):
        return [self.get_model().get_particle(pid) for pid in self.pids]



# Return a list of restraints rs on the IMP model m (with hierarchy h) based on potential energy terms from the CHARMM force-field/parameterization.
def charmm_restraints(
        m,
        h,
        eps=False
):

    rs = list()

    # Configure the IMP model based on the CHARMM parameterization.
    # ff = IMP.atom.get_heavy_atom_CHARMM_parameters()
    # pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    # print(len(pids))

    ff = IMP.atom.get_all_atom_CHARMM_parameters()
    topology = ff.create_topology(h)

    # topology.apply_default_patches()
    # topology.setup_hierarchy(h)
    topology.add_atom_types(h)
    # topology.add_missing_atoms(h)
    # IMP.atom.remove_charmm_untyped_atoms(h)
    # topology.add_coordinates(h)

    # pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    # print(len(pids))


    # topology.add_missing_atoms(h)
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