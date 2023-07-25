from pathlib import Path

import IMP
import IMP.atom


def charmmm_met(
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

    improper_p = impropers[0]
    print(m.get_particle_name(improper_p))
    print(IMP.atom.Dihedral.get_is_setup(m, improper_p.get_index()))

    d_p = IMP.atom.Dihedral(m, improper_p.get_index())
    d_0_p = d_p.get_particle(0)

    print(IMP.atom.Atom.get_is_setup(m, d_0_p.get_index()))
    print(m.get_particle_name(d_0_p.get_index()))

    pids_235 = IMP.atom.Selection(hierarchy=h, residue_index=235)
    pids_main_chain = IMP.atom.Selection(hierarchy=h, atom_types=[IMP.atom.AT_CA, IMP.atom.AT_C, IMP.atom.AT_O, IMP.atom.AT_N])

    pids_235_side_chain = pids_235-pids_main_chain

    containers = {"bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": impropers}
    containers_sub = dict()
    for container_name in containers.keys():
        print(container_name, len(containers[container_name]))

    for container_name in containers.keys():
        print(container_name)

        # The general strategy for removing bonds is by iterating through all of the atoms in the side chains list and then deleting any bond containing that atom.
        if container_name == "bonds":
            bonds_copy = bonds.copy()
            bond_indexes = list()

            for atom_pid in pids_235_side_chain.get_selected_particle_indexes():
                bonded = IMP.atom.Bonded(m, atom_pid)
                print(m.get_particle_name(atom_pid))
                print(bonded.get_number_of_bonds())

                bonded_bonds = [bonded.get_bond(i) for i in range(bonded.get_number_of_bonds())]

                for bond in bonded_bonds:
                    print(bond in bonds, bond in bonds_copy)
                    bond_index = 0
                    for p in bonds_copy:
                        if bond == p:
                            break
                        bond_index = bond_index + 1

                    print(bond_index)
                    bond_indexes.append(bond_index)

            bond_indexes = sorted(list(set(bond_indexes)))
            print(bond_indexes)

            # Deleting elements
            offset = 0
            for bond_index in bond_indexes:
                del bonds_copy[bond_index-offset]
                offset = offset+1

            containers_sub["bonds"] = bonds_copy

    for container_name in containers_sub.keys():
        print(container_name, len(containers_sub[container_name]))


        # ps = containers[container_name]
        # for p in ps:
        #     pid = p.get_index()

        #     if container_name in ["bonds"]:
                # bond_dec = IMP.atom.Bond(m, pid)

                # # bonded_1 and bonded_2 are decorators of type IMP.atom.bonded.
                # bonded_1 = bond_dec.get_bonded(0)
                # bonded_2 = bond_dec.get_bonded(1)

                # atom_pids = []




            # elif container_name in ["angles"]:
            #     decorator = IMP.atom.Angle(m, pid)
            #     n_atoms = 3

            # else:
            #     decorator = IMP.atom.Dihedral(m, pid)
            #     atom_pids = [decorator.get_particle(i).get_index() for i in range(n_atoms)]

            #     n_atoms = 4


            # contains = False
            # for atom_pid in atom_pids:
            #     if atom_pid in pids_235_side_chain.get_selected_particle_indexes():
            #         contains = True

            # if contains:
            #     print(m.get_particle_name(pid))
            #     for atom_pid in atom_pids:
            #         print(m.get_particle_name(atom_pid))

        break

    # for dihedral in dihedrals:
    #     pid = dihedral.get_index()
    #     if pid in pids_235_side_chain.get_selected_particle_indexes():
    #         print(m.get_particle_name(dihedral))

    # for bond in bonds:
    #     pid = bond.get_index()
    #     if pid in pids_235_side_chain.get_selected_particle_indexes():
    #         print(m.get_particle_name(bond))

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

    return rs


if __name__ == "__main__":
    pdb_file = Path("/wynton/group/sali/mhancock/xray/dev/12_sulfur_only/data/out/01_300/9368668/output_4/pdbs/100.pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    rs = list()
    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmmm_met(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)
    rset_charmm.set_weight(1)
    rs.append(rset_charmm)