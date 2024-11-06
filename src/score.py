import sys
from pathlib import Path
import multiprocessing
import itertools
import pandas as pd
import numpy as np

import IMP
import IMP.atom
import IMP.isd

from align_imp import get_multi_state_multi_cond_rmsd
import cctbx_score
import charmm
import miller_ops
import params
import weights
from multi_state_multi_condition_model import MultiStateMultiConditionModel
from xray_restraint import XtalRestraint
from utility import get_n_state_from_pdb_file
from miller_ops import get_f_obs, get_flags
from charmm import get_charmm_restraint_set


"""
The function "score" takes in various parameters such as hs, hs_0, cif_file, flags_file, res, and score_fs, and calculates different scores based on the given parameters. It uses the CHARMM program to compute force field (ff) scores and the CCTBX library to calculate maximum likelihood (ml) and least squares (ls) scores. The function also computes the root-mean-square deviation (rmsd) between two sets of input data. The scores are stored in a dictionary with the corresponding score term as the key.

Parameters:
    hs: A list of hierarchies. The occupancies of these heirachies will be used in score calculations.

    hs_0: A list of reference hierachies used for computing RMSD.

    f_obs: A miller array containing the observed structure factor amplitudes.

    flags: A miller array containing the r_free flags for the observed structure factor amplitudes.

    score_fs: A list of score terms for which scores need to be calculated.

Return:
    scores_dict: A dictionary containing the calculated scores, with score terms as keys.
"""
def pool_score(
        params
):
    decoy_file=params["decoy_file"]
    decoy_w_mat = params["decoy_w_mat"]
    ref_file=params["ref_file"]
    ref_w_mat=params["ref_w_mat"]
    score_fs=params["score_fs"]

    if "cif_files" in params:
        cif_files=params["cif_files"]
    else:
        cif_files=None

    scores_dict = dict()
    # Indicate the native structure.
    if decoy_file == ref_file:
        scores_dict["native"] = 1
    else:
        scores_dict["native"] = 0

    scores_dict["pdb_file"] = str(decoy_file)

    # decoy_w_mat = np.ndarray(shape=[len(decoy_occs), 1])
    # decoy_w_mat[:,0] = decoy_occs

    if cif_files:
        crystal_symmetries = list()
        for cif_file in cif_files:
            crystal_symmetries.append(miller_ops.get_crystal_symmetry(cif_file))
    else:
        crystal_symmetries = None

    decoy_msmc_m = MultiStateMultiConditionModel(
        pdb_files=[decoy_file],
        w_mat=decoy_w_mat,
        crystal_symmetries=crystal_symmetries
    )
    h_decoys = decoy_msmc_m.get_hs()

    # ref_w_mat = np.ndarray(shape=[len(ref_occs), 1])
    # ref_w_mat[:,0] = ref_occs
    ref_msmc_m = MultiStateMultiConditionModel(
        pdb_files=[ref_file],
        w_mat=ref_w_mat,
        crystal_symmetries=None
    )

    for score_f in score_fs:
        if score_f in ["ff", "bnd", "ang", "dih", "imp", "eps", "nbd"]:
            charmm_rset = get_charmm_restraint_set(decoy_msmc_m.get_m(), h_decoys)
            score = charmm_rset.evaluate(False)
        elif "xray" in score_f:
            cond = int(score_f.split("_")[1])
            cif_file = cif_files[cond]

            f_obs = get_f_obs(cif_file)
            flags = get_flags(cif_file)
            f_obs, flags = f_obs.common_sets(other=flags)

            res = params["res"]
            # print(f_obs.size())
            if res:
                f_obs = miller_ops.filter_f_obs_resolution(
                    f_obs=f_obs,
                    d_max=None,
                    d_min=res
                )
                flags = miller_ops.filter_f_obs_resolution(
                    f_obs=flags,
                    d_max=None,
                    d_min=res
                )
            # print(f_obs.size())

            pids = list()
            for h in h_decoys:
                pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

            xray_r = XtalRestraint(
                msmc_m=decoy_msmc_m,
                cond=cond,
                f_obs=f_obs,
                free_flags=flags,
                w_xray=1,
                update_scale=params["scale"],
                update_k1=params["scale_k1"],
                ref_com=ref_msmc_m.get_com(),
                update_freq=1
            )
            xray_r.evaluate(False)
            score = xray_r.get_f()
            scores_dict["r_free_{}".format(cond)] = xray_r.get_r_free()
            scores_dict["r_work_{}".format(cond)] = xray_r.get_r_work()

            scores_dict["cif_files"] = cif_files
        elif "rmsd" == score_f:
            score = get_multi_state_multi_cond_rmsd(
                decoy_msmc_m,
                ref_msmc_m
            )
        elif "rmsd" in score_f:
            cond = int(score_f.split("_")[1])
            score = get_multi_state_multi_cond_rmsd(
                decoy_msmc_m,
                ref_msmc_m,
                cond
            )

            # if score_f == "rmsd_avg":
            #     f = align_imp.compute_rmsd_between_average
            # elif score_f == "rmsd_ord":
            #     f = align_imp.compute_rmsd_ordered
            # elif score_f == "rmsd_dom":
            #     f = align_imp.compute_rmsd_dom_state
            # elif score_f == "avg_delta_w":
            #     f = align_imp.compute_avg_delta_weight

            # ## NEED TO FIX
            # score = f(
            #         h_0s=h_decoys,
            #         h_1s=ref_msmc_m.get_hs(),
            #         pids_0=decoy_msmc_m.get_ca_pids(),
            #         pids_1=ref_msmc_m.get_ca_pids(),
            #         occs_0=decoy_msmc_m.get_occs_for_condition_i(0),
            #         occs_1=ref_msmc_m.get_occs_for_condition_i(0),
            #     )
        elif score_f == "w":
            score = list(w.get_weights())
        else:
            score = np.nan

        scores_dict[score_f] = score

    return scores_dict

