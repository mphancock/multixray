
def add_H_to_N_ter(
    m, h,
    topology
):
    res_id = 1
    pid = IMP.atom.Selection(h,residue_index=res_id,atom_type=IMP.atom.AtomType("H1")).get_selected_particle_indexes()[0]
    IMP.atom.CHARMMAtom.setup_particle(m, pid, "HC")
    charge = IMP.atom.Charged.setup_particle(m, pid, 0.33)

    at = IMP.atom.Atom(m, pid)
    charmm_at = IMP.atom.CHARMMAtom(m, pid)
    print("PATCHING: ", at.get_name(), charmm_at.get_charmm_type())


    segment = topology.get_segment(0)
    res = segment.get_residue(res_id-1)
    res.set_patched(False)

    ## NEED TO USE PDB NAME NOT CHARMM NAME
    bond = IMP.atom.CHARMMBond(["H1", "N"])
    patch = IMP.atom.CHARMMPatch("TMP")
    patch.add_bond(bond)
    patch.apply(res)


## need to change histidine from HIS to HSP
def convert_HIS_to_HSP(
    m, h,
    ff,
    topology,
    res_id
):
    seg = topology.get_segments()[0]
    ress = seg.get_residues()
    # for res_id in [41, 80, 163, 172, 246]:

    ## check if the histidine contains an HE2 atom
    pids = IMP.atom.Selection(h,residue_index=res_id,atom_type=IMP.atom.AtomType("HE2")).get_selected_particle_indexes()

    if len(pids) == 0:
        return

    print("PATCHING HIS: ", res_id)

    he2_pid = pids[0]
    ## setup the HE2 (doesn't exist in CHARMM) atom as H
    IMP.atom.CHARMMAtom.setup_particle(m, he2_pid, "H")
    charge = IMP.atom.Charged.setup_particle(m, he2_pid, 0.09)

    ## update the topology
    ress[res_id-1] = IMP.atom.CHARMMResidueTopology(ff.get_residue_topology(IMP.atom.ResidueType("HSP")))

    ## change the charmm type of the N from NR2 (unprotonated) to NR3 (protonated)
    ne2_pid = IMP.atom.Selection(h, residue_index=res_id, atom_type=IMP.atom.AtomType("NE2")).get_selected_particles()[0]
    IMP.atom.CHARMMAtom(m, ne2_pid).set_charmm_type("NR3")

    seg.set_residues(ress)
