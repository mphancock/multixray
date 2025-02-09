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
    h_0s (list): list of hierarchies for the first model.
    h_1s (list): list of hierarchies for the second model.
    pids_0 (list): list of particle ids for one state of the first model.
    pids_1 (list): list of particle ids for one state of the second model.
    occs_0 (list): list of weights for the first model.
    occs_1 (list): list of weights for the second model.

Returns
    rmsd (float): the RMSD between the two models.
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
Function to compute the RMSDs between the states of 2 multi-state multi-condition models. The states are ordered based on whatever condition is being compared. The 2 multi-state multi-condition models must contain the same number of compositionally identical states.

Params
    msmc_0 (MultiStateMultiConditionModel): the first multi-state multi-condition model.
    msmc_1 (MultiStateMultiConditionModel): the second multi-state multi-condition model.
    cond (int): the condition to compare.

Returns
    rmsds (list): the RMSDs between the states of the two models.
"""
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


"""
Function to compute the weight errors between the states of 2 multi-state multi-condition models. The states are ordered based on whatever condition is being compared. The 2 multi-state multi-condition models must contain the same number of compositionally identical states.

Params
    msmc_0 (MultiStateMultiConditionModel): the first multi-state multi-condition model.
    msmc_1 (MultiStateMultiConditionModel): the second multi-state multi-condition model.
    cond (int): the condition to compare.

Returns
    errors (list): the weight errors between the states of the two models.
"""
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


def align_one_to_two(
    h_1,
    h_2
):
    m_1 = h_1.get_model()
    m_2 = h_2.get_model()

    ca_pids_1 = IMP.atom.Selection(h_1, atom_type=IMP.atom.AtomType("CA")).get_selected_particles()
    ca_pids_2 = IMP.atom.Selection(h_2, atom_type=IMP.atom.AtomType("CA")).get_selected_particles()

    xyzs_1 = [IMP.core.XYZ(m_1, ca_pid).get_coordinates() for ca_pid in ca_pids_1]
    xyzs_2 = [IMP.core.XYZ(m_2, ca_pid).get_coordinates() for ca_pid in ca_pids_2]

    transformation = IMP.algebra.get_transformation_aligning_first_to_second(xyzs_1, xyzs_2)
    transform = IMP.core.Transform(transformation)
    transform.apply_indexes(m_1, IMP.atom.Selection(h_1).get_selected_particle_indexes(), 0, 9999)