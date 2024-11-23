from pathlib import Path
import numpy as np

import IMP
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.isd

import average_structure
from utility import get_ordered_hs


"""
Function to compute the RMSD between 2 multi-state models. The function assumes that the 2 models contain compositionally identical states. The multi-state RMSD is computed by averaging the coordinates of the atoms in each state and then computing the RMSD between the two average structures.

Params
    h_0s: list of hierarchies for the first model.
    h_1s: list of hierarchies for the second model.
    pids_0: list of particle ids for one state of the first model.
    pids_1: list of particle ids for one state of the second model.
    occs_0: list of weights for the first model.
    occs_1: list of weights for the second model.

Returns
    rmsd: the RMSD between the two models.
"""
def get_multi_state_rmsd(
    h_0s,
    h_1s,
    pids_0,
    pids_1,
    occs_0,
    occs_1
):
    avg_dict_0 = average_structure.get_coord_avg_dict(
        hs=h_0s,
        occs=occs_0
    )

    avg_dict_1 = average_structure.get_coord_avg_dict(
        hs=h_1s,
        occs=occs_1
    )

    ## check that the length of the dictionaries is equal
    if len(avg_dict_0) != len(avg_dict_1):
        raise RuntimeError("The number of atoms in the average structures are not equal: {} and {}".format(len(avg_dict_0), len(avg_dict_1)))

    rmsd = 0
    for i in range(len(pids_0)):
        pid_0 = pids_0[i]
        pid_1 = pids_1[i]

        name_0 = h_0s[0].get_model().get_particle_name(pid_0)
        name_1 = h_1s[0].get_model().get_particle_name(pid_1)

        if name_0 != name_1:
            raise RuntimeError("Particle names not equal: {} and {}".format(name_0, name_1))

        xyz_1 = avg_dict_0[pids_0[i]]
        xyz_2 = avg_dict_1[pids_1[i]]
        mag = (xyz_1-xyz_2).get_magnitude()

        rmsd = rmsd+mag**2

    rmsd = rmsd/len(pids_0)
    rmsd = np.sqrt(rmsd)

    return float(rmsd)


"""
Function to compute the RMSD between 2 multi-state multi-condition models. The function assumes that the 2 models contain compositionally identical states and equal number of conditions.

Params
    msmc_0: the first multi-state multi-condition model.
    msmc_1: the second multi-state multi-condition model.

Returns
    rmsd: the RMSD between the two models.
"""
def get_multi_state_multi_cond_rmsd(
    msmc_0,
    msmc_1,
    cond=None # None means all conditions
):
    if msmc_0.get_n_cond() != msmc_1.get_n_cond():
        raise RuntimeError("Number of conditions not equal: {} and {}".format(msmc_0.get_n_cond(), msmc_1.get_n_cond()))

    rmsd = 0
    # for cond in range(msmc_0.get_n_cond()):

    if cond is None:
        conds = range(msmc_0.get_n_cond())
    else:
        conds = [cond]

    for cond in conds:
        cond_rmsd = get_multi_state_rmsd(
            h_0s=msmc_0.get_hs(),
            h_1s=msmc_1.get_hs(),
            pids_0=msmc_0.get_ca_pids_in_state(0),
            pids_1=msmc_1.get_ca_pids_in_state(0),
            occs_0=msmc_0.get_occs_for_condition(cond).tolist(),
            occs_1=msmc_1.get_occs_for_condition(cond).tolist()
        )
        rmsd += cond_rmsd

    return rmsd / len(conds)


"""
Function to compute the RMSD between 2 pdb files containing multi-state models. The function assumes that the 2 models contain compositionally identical states.

Params
    pdb_0: the first pdb file.
    pdb_1: the second pdb file.

Returns
    rmsd: the RMSD between the two models.
"""
def get_multi_state_rmsd_from_pdbs(
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


def compute_rmsds_between_ordered_states(
        msmc_0,
        msmc_1,
        cond
):
    if msmc_0.get_n_state() != msmc_1.get_n_state():
        raise RuntimeError("Number of states not equal: {} and {}".format(len(h_0s), len(h_1s)))

    h_0s = get_ordered_hs(msmc_0.get_hs(), msmc_0.get_occs_for_condition(cond))
    h_1s = get_ordered_hs(msmc_1.get_hs(), msmc_1.get_occs_for_condition(cond))

    rmsds = list()
    for i in range(len(h_0s)):
        sel_0 = IMP.atom.Selection(h_0s[i], atom_type=IMP.atom.AtomType("CA"))
        sel_1 = IMP.atom.Selection(h_1s[i], atom_type=IMP.atom.AtomType("CA"))
        rmsds.append(IMP.atom.get_rmsd(sel_0, sel_1))

    return rmsds


def compute_weight_errors_between_ordered_states(
        msmc_0,
        msmc_1,
        cond
):
    if msmc_0.get_n_state() != msmc_1.get_n_state():
        raise RuntimeError("Number of states not equal: {} and {}".format(len(h_0s), len(h_1s)))

    h_0s = get_ordered_hs(msmc_0.get_hs(), msmc_0.get_occs_for_condition(cond))
    h_1s = get_ordered_hs(msmc_1.get_hs(), msmc_1.get_occs_for_condition(cond))

    occs_0 = msmc_0.get_occs_for_condition(cond).tolist()
    occs_0.sort(reverse=True)
    occs_1 = msmc_1.get_occs_for_condition(cond).tolist()
    occs_1.sort(reverse=True)

    errors = list()
    for i in range(len(h_0s)):
        occ_0 = occs_0[i]
        occ_1 = occs_1[i]

        errors.append((occ_0-occ_1)**2)

    return errors


if __name__ == "__main__":
    from multi_state_multi_condition_model import MultiStateMultiConditionModel

    decoy_msmc_m = MultiStateMultiConditionModel(
        pdb_files=[Path("/wynton/home/sali/mhancock/xray/dev/45_synthetic_native_4/data/pdbs/native.pdb")],
        w_mat=np.array([[1, 0], [0, 1]]),
        crystal_symmetries=None
    )
    h_decoys = decoy_msmc_m.get_hs()

    # ref_w_mat = np.ndarray(shape=[len(ref_occs), 1])
    # ref_w_mat[:,0] = ref_occs
    # ref_msmc_m = MultiStateMultiConditionModel(
    #     pdb_files=[Path("/wynton/home/sali/mhancock/xray/dev/45_synthetic_native_4/data/pdbs/native.pdb")],
    #     w_mat=np.array([[.9, .1], [.1, 0.9]]),
    #     crystal_symmetries=None
    # )

    dict_0 = average_structure.get_coord_avg_dict(
        hs=decoy_msmc_m.get_hs(),
        occs=[1, 0]
    )

    dict_1 = average_structure.get_coord_avg_dict(
        hs=decoy_msmc_m.get_hs(),
        occs=[0, 1]
    )

    test_pid = decoy_msmc_m.get_ca_pids_in_state(0)[0]
    print(dict_0[test_pid])
    print(dict_1[test_pid])

    # print(get_multi_state_multi_cond_rmsd(decoy_msmc_m, ref_msmc_m, 0))
    # print(get_multi_state_multi_cond_rmsd(decoy_msmc_m, ref_msmc_m, 1))