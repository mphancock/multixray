import sys
from pathlib import Path

import IMP
import IMP.atom
import IMP.core
import IMP.algebra

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))


if __name__ == "__main__":
    m = IMP.Model()
    pdb_file = Path("../data/min.pdb")

    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    params_file = Path("/Users/matthew/opt/anaconda3/envs/imp_221_cctbx/share/IMP/atom/par.lib")
    topology_file = Path("/Users/matthew/opt/anaconda3/envs/imp_221_cctbx/share/IMP/atom/top.lib")

    # ff = IMP.atom.CHARMMParameters(str(topology_file), str(params_file), True)
    # topology = ff.create_topology(h)

    ff = IMP.atom.get_all_atom_CHARMM_parameters()
    topology = ff.create_topology(h)

    topology.apply_default_patches()
    topology.add_atom_types(h)


    rs = list()

    # Configure the IMP model based on the CHARMM parameterization.
    # ff = IMP.atom.get_heavy_atom_CHARMM_parameters()
    # pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    # print(len(pids))

    ff = IMP.atom.get_all_atom_CHARMM_parameters()
    topology = ff.create_topology(h)

    ## Why is this not fixing the NTER?
    topology.apply_default_patches()
    # topology.setup_hierarchy(h)
    topology.add_atom_types(h)
    # topology.add_missing_atoms(h)
    IMP.atom.remove_charmm_untyped_atoms(h)
    # topology.add_coordinates(h)

    ## setup charges first
    charges = topology.add_charges(h)
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

    slack = 0

    # Add a restraint on the non-bonded atoms (Lennard-Jones potential).
    ff.add_radii(h)
    ff.add_well_depths(h)
    atoms = IMP.atom.get_by_type(h, IMP.atom.ATOM_TYPE)
    cont = IMP.container.ListSingletonContainer(m, atoms)
    nbl = IMP.container.ClosePairContainer(cont, 3, slack)
    pair_filter = IMP.atom.StereochemistryPairFilter()
    pair_filter.set_bonds(bonds)
    pair_filter.set_angles(angles)
    pair_filter.set_dihedrals(dihedrals)
    nbl.add_pair_filter(pair_filter)
    sf = IMP.atom.ForceSwitch(6.0, 7.0)
    ljps = IMP.atom.LennardJonesPairScore(sf)
    rs.append(IMP.container.PairsRestraint(ljps, nbl, "nbd"))

    sf = IMP.atom.ForceSwitch(6.0, 7.0)
    cps = IMP.atom.CoulombPairScore(sf)
    r_cps = IMP.container.PairsRestraint(cps, nbl, "eps")
    rs.append(r_cps)

    pid = IMP.atom.Selection(h).get_selected_particle_indexes()[0]

    for r in rs:
        print(r.get_name(), r.evaluate(True))
        print(IMP.core.XYZ(m, pid).get_derivatives())

    print(IMP.core.XYZ(m, pid).get_derivatives())

    pids = IMP.atom.Selection(h).get_selected_particles()
    ps = [m.get_particle(pid) for pid in pids]
    for pid in pids:
        IMP.core.XYZ(m, pid).set_coordinates_are_optimized(True)
        IMP.atom.LinearVelocity.setup_particle(m, pid, IMP.algebra.Vector3D(0, 0, 0))

    sf = IMP.core.RestraintsScoringFunction(rs)

    sf.evaluate(True)

    pid = IMP.atom.Selection(h).get_selected_particle_indexes()[0]

    # vel_therm = IMP.atom.BerendsenThermostatOptimizerState(ps, 300, 10)

    md = IMP.atom.MolecularDynamics(m)

    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    md.setup(ps)
    md.set_temperature(300)
    md.set_maximum_time_step(1.0)
    md.assign_velocities(300)

    # md.add_optimizer_state(vel_therm)

    md.simulate(1000)


    for r in rs:
        print(r.get_name(), r.evaluate(True))
        print(IMP.core.XYZ(m, pid).get_derivatives())

    print(IMP.core.XYZ(m, pid).get_derivatives())