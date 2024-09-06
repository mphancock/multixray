from pathlib import Path
import numpy as np

import IMP
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.isd

import average_structure



"""
Given 2 sets of IMP hierarchies, the function returns the RMSD of the average structure of each set of hierarchies. This allows RMSD to be computed between sets of hierarchies of unequal size (eg, a decoy multi-state structure against a single-state native). One challenge is that is important that any multi state structure contain identical copies otherwise average structures cannot be computed.

Consider the following situation. We have a 1-state set of structures (h1) and a 2-state set of structures (h21, h22). If you try and compute the average structure between [h1] and [h21, h22], it will succeed becasue the pids of the average structures will match (because h1 and h21 match). If you try and compute an average structure between [h1] and [h22] it will fail because the average [h22] structure will have different set of pids. In this situation, we cannot rely on pids.

**********
Params
    h_1s: list of IMP hierarchies of size N.

    h_2s: list of IMP hierarchies of size M.

    ca_only: boolean, True if you only want to compute the RMSD between CA atoms of the hierarchies.

**********
Returns
    rmsd: the rmsd between the average structures of h_1s and h_2s.

"""
def compute_rmsd_between_average(
        h_0s,
        h_1s,
        pids_0,
        pids_1,
        occs_0,
        occs_1
):
    avg_dict_1 = average_structure.get_coord_avg_dict(
        hs=h_0s,
        occs=occs_0
    )

    avg_dict_2 = average_structure.get_coord_avg_dict(
        hs=h_1s,
        occs=occs_1
    )

    rmsd = 0
    for i in range(len(pids_0)):
        pid_0 = pids_0[i]
        pid_1 = pids_1[i]

        name_0 = h_0s[0].get_model().get_particle_name(pid_0)
        name_1 = h_1s[0].get_model().get_particle_name(pid_1)

        if name_0 != name_1:
            raise RuntimeError("Particle names not equal: {} and {}".format(name_0, name_1))

        # print(h_0s[0].get_model().get_particle_name(pids_0[i]), h_1s[0].get_model().get_particle_name(pids_1[i]))
        xyz_1 = avg_dict_1[pids_0[i]]
        xyz_2 = avg_dict_2[pids_1[i]]
        mag = (xyz_1-xyz_2).get_magnitude()
        rmsd = rmsd+mag**2

    rmsd = rmsd/len(pids_0)
    rmsd = np.sqrt(rmsd)

    return rmsd


def compute_rmsd_between_average_pdb(
        pdb_0,
        pdb_1
):
    m, m_0 = IMP.Model(), IMP.Model()
    hs_0 = IMP.atom.read_multimodel_pdb(str(pdb_0), m, IMP.atom.AllPDBSelector())
    hs_1 = IMP.atom.read_multimodel_pdb(str(pdb_1), m_0, IMP.atom.AllPDBSelector())

    pids_0 = IMP.atom.Selection(hs_0[0]).get_selected_particle_indexes()
    pids_1 = IMP.atom.Selection(hs_1[0]).get_selected_particle_indexes()

    occs_0 = np.array([1/len(hs_0)]*len(hs_0))
    occs_1 = np.array([1/len(hs_1)]*len(hs_1))

    rmsd = compute_rmsd_between_average(
        h_0s=hs_0,
        h_1s=hs_1,
        pids_0=pids_0,
        pids_1=pids_1,
        occs_0=occs_0,
        occs_1=occs_1
    )

    return rmsd


def compute_weighted_rmsd(
        h_0s,
        h_1s,
        ca_only
):
    n_state = len(h_0s)
    if len(h_0s) != len(h_1s):
        raise RuntimeError("List of hierarchies not of equal sizes: {} and {}".format(len(h_0s), len(h_1s)))

    if ca_only:
        pids = IMP.atom.Selection(h_0s[0], atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()
    else:
        pids = IMP.atom.Selection(h_0s[0]).get_selected_particle_indexes()

    rmsd = 0
    for i in range(len(h_0s)):
        h_1 = h_0s[i]
        h_2 = h_1s[i]

        for pid in pids:
            coord_1 = IMP.core.XYZ(h_1.get_model(), pid).get_coordinates()
            coord_2 = IMP.core.XYZ(h_2.get_model(), pid).get_coordinates()
            occ_1 = IMP.atom.Atom(h_1.get_model(), pid).get_occupancy()
            occ_2 = IMP.atom.Atom(h_2.get_model(), pid).get_occupancy()

            mag = (coord_1*occ_1-coord_2*occ_2).get_magnitude()
            rmsd = rmsd+mag**2

    rmsd = rmsd/(len(pids)*n_state)
    rmsd = np.sqrt(rmsd)

    return rmsd


# Compute the root mean squared deviation between 2 IMP models.
def compute_rmsd(
        h_0s,
        h_1s,
        ca_only
):
    n_state = len(h_0s)
    rmsd = 0
    for i in range(n_state):
        h_0 = h_0s[i]
        h_1 = h_1s[i]

        if ca_only:
            pids_0 = IMP.atom.Selection(h_0, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()
            pids_1 = IMP.atom.Selection(h_1, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()
        else:
            pids_0 = IMP.atom.Selection(h_0).get_selected_particle_indexes()
            pids_1 = IMP.atom.Selection(h_1).get_selected_particle_indexes()

        for j in range(len(pids_0)):
            pid_0 = pids_0[j]
            pid_1 = pids_1[j]

            coord_0 = IMP.core.XYZ(h_0.get_model(), pid_0).get_coordinates()
            coord_1 = IMP.core.XYZ(h_1.get_model(), pid_1).get_coordinates()

            occ_0 = IMP.atom.Atom(h_0.get_model(), pid_0).get_occupancy()

            mag = (coord_0-coord_1).get_magnitude()
            rmsd = rmsd + occ_0 * mag**2

    rmsd = rmsd/(len(pids_0))
    rmsd = np.sqrt(rmsd)

    return rmsd


def get_ordered_hs(
        h_0s
):
    occs = list()
    for i in range(len(h_0s)):
        h_0 = h_0s[i]
        pid_0 = IMP.atom.Selection(h_0).get_selected_particles()[0]
        occ = IMP.atom.Atom(h_0.get_model(), pid_0).get_occupancy()
        occs.append(occ)

    # Need to invert the list.
    ids = list(np.argsort(occs))
    ids.reverse()

    h_0s_ordered = [h_0s[i] for i in ids]

    return h_0s_ordered


def set_occupancies(
    h,
    occ
):
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    for pid in pids:
        IMP.atom.Atom(h.get_model(), pid).set_occupancy(occ)


def single_h_to_hs(
    h,
    n_state
):
    h_clones = list()

    for i in range(n_state):
        h_clone = IMP.atom.create_clone(h)
        h_clones.append(h_clone)

    for h in h_clones:
        set_occupancies(h=h, occ=1/n_state)

    return h_clones


def compute_rmsd_ordered(
        h_0s,
        h_1s,
        ca_only
):
    clone_0, clone_1 = False, False

    if len(h_0s) != len(h_1s):
        n_state = max(len(h_0s), len(h_1s))
        if len(h_0s) == 1:
            h_0s = single_h_to_hs(h=h_0s[0], n_state=n_state)
            clone_0 = True
        elif len(h_1s) == 1:
            h_1s = single_h_to_hs(h=h_1s[0], n_state=n_state)
            clone_1 = True
        else:
            raise RuntimeError("h_0s and h_1s not of equal size: {} and {}".format(len(h_0s), len(h_1s)))

    h_0s_ordered = get_ordered_hs(h_0s)
    h_1s_ordered = get_ordered_hs(h_1s)

    rmsd = compute_multi_rmsd(
        h_0s=h_0s_ordered,
        h_1s=h_1s_ordered,
        ca_only=ca_only
    )

    # Need to cleanup the hierarchies because they are attatched to the hierarchies underlying model.
    if clone_0:
        for h in h_0s:
            IMP.atom.destroy(h)
    elif clone_1:
        for h in h_1s:
            IMP.atom.destroy(h)

    return rmsd


def compute_multi_rmsd(
        h_0s,
        h_1s,
        ca_only
):
    if len(h_0s) != len(h_1s):
        raise RuntimeError("h_0s and h_1s not of equal size: {} and {}".format(len(h_0s), len(h_1s)))

    n_state = len(h_0s)
    rmsd_tot = 0

    for i in range(n_state):
        if ca_only:
            rmsd = IMP.atom.get_rmsd(IMP.atom.Selection(h_0s[i], atom_type=IMP.atom.AtomType("CA")), IMP.atom.Selection(h_1s[i], atom_type=IMP.atom.AtomType("CA")))
        else:
            rmsd = IMP.atom.get_rmsd(IMP.atom.Selection(h_0s[i]), IMP.atom.Selection(h_1s[i]))

        occ_0 = IMP.atom.Atom(h_0s[i].get_model(), IMP.atom.Selection(h_0s[i]).get_selected_particle_indexes()[0]).get_occupancy()

        rmsd_tot = rmsd_tot + occ_0 * rmsd

    return rmsd


def compute_rmsd_dom_state(
        h_0s,
        h_1s,
        ca_only
):
    h_0 = get_ordered_hs(h_0s)[0]
    h_1 = get_ordered_hs(h_1s)[0]

    if ca_only:
        rmsd = IMP.atom.get_rmsd(IMP.atom.Selection(h_0, atom_type=IMP.atom.AtomType("CA")), IMP.atom.Selection(h_1, atom_type=IMP.atom.AtomType("CA")))
    else:
        rmsd = IMP.atom.get_rmsd(IMP.atom.Selection(h_0), IMP.atom.Selection(h_1))

    return rmsd


def compute_avg_delta_weight(
        h_0s,
        h_1s,
        ca_only
):
    n_states = len(h_0s)
    delta = 0

    h_0s_ordered = get_ordered_hs(h_0s)
    h_1s_ordered = get_ordered_hs(h_1s)

    for i in range(n_states):
        h_0 = h_0s_ordered[i]
        h_1 = h_1s_ordered[i]

        pid_0 = IMP.atom.Selection(h_0).get_selected_particles()[0]
        pid_1 = IMP.atom.Selection(h_1).get_selected_particles()[0]

        occ_0 = IMP.atom.Atom(h_0.get_model(), pid_0).get_occupancy()
        occ_1 = IMP.atom.Atom(h_1.get_model(), pid_1).get_occupancy()

        delta = delta + np.sqrt((occ_0-occ_1)**2)

    avg_delta = delta / n_states

    return avg_delta


if __name__ == "__main__":
    pdb_selector = IMP.atom.NonAlternativePDBSelector()

    pdb_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30/0.pdb")
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, pdb_selector)
    pids = IMP.atom.Selection(hs[0], atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()
    occs = np.array([0.510587948937655,0.489412051062345])

    ref_pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/121_native_decoys/0/output_24/pdbs/7.pdb")
    ref_m = IMP.Model()
    ref_hs = IMP.atom.read_multimodel_pdb(str(ref_pdb_file), ref_m, IMP.atom.AllPDBSelector())
    ref_pids = IMP.atom.Selection(ref_hs[0], atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()
    ref_occs = np.array([0.4098692178088230, 0.5901307821911770])

    # ref_pdb_file = Path("/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhk.pdb")
    # ref_m = IMP.Model()
    # ref_hs = IMP.atom.read_multimodel_pdb(str(ref_pdb_file), ref_m, pdb_selector)
    # ref_pids = IMP.atom.Selection(ref_hs[0], atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()
    # ref_occs = np.array([1])

    # rmsd = compute_rmsd_between_average(
    #     h_0s=hs,
    #     h_1s=ref_hs,
    #     pids_0=pids,
    #     pids_1=ref_pids,
    #     occs_0=occs,
    #     occs_1=ref_occs
    # )

    print(len(pids), len(ref_pids))
    print(rmsd)