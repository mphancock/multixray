import sys
from pathlib import Path
import IMP
import IMP.atom

import IMP
import IMP.core
import IMP.algebra


if __name__ == "__main__":
    ## IMP.get_data_path("charmm/top_all22_prot.inp")
    m = IMP.Model()
    pdb_file = Path(Path.home(), "Documents/xray/data/pdbs/7mhf/7mhf.pdb")

    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    # params_file = Path("/Users/matthew/opt/anaconda3/envs/imp_221_cctbx/share/IMP/atom/par.lib")
    # topology_file = Path("/Users/matthew/opt/anaconda3/envs/imp_221_cctbx/share/IMP/atom/top.lib")
    params_file = Path("/Users/matthew/Documents/xray/tmp/params.txt")
    topology_file = Path("/Users/matthew/Documents/xray/tmp/topology.txt")

    ff = IMP.atom.CHARMMParameters(str(topology_file), str(params_file), True)

    print(ff.get_bond_parameters("NH3", "HC"))

    topology = ff.create_topology(h)

    topology.apply_default_patches()
    topology.add_atom_types(h)
    # topology.add_missing_atoms(h)

    rs = list()
    # topology.add_missing_atoms(h)
    bonds = topology.add_bonds(h)
    angles = ff.create_angles(bonds)
    dihedrals = ff.create_dihedrals(bonds)
    impropers = topology.add_impropers(h)
    # charges = topology.add_charges(h)

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

    sf = IMP.core.RestraintsScoringFunction(rs)
    cg = IMP.core.ConjugateGradients(m)
    cg.set_scoring_function(sf)

    cg.optimize(100)

    out_file = Path(Path.home(), "Documents/xray/tmp/tmp.pdb")
    IMP.atom.write_pdb(h, str(out_file))