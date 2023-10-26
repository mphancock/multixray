import sys
from pathlib import Path
import multiprocessing
import itertools
import pandas as pd
import numpy as np

import IMP
import IMP.atom
import IMP.isd

sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp
# import cctbx_scores
import cctbx_score
import charmm
import miller_ops
import params


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
    occs=params["occs"]
    ref_file=params["ref_file"]
    adp_file=params["adp_file"]
    cif_file=params["cif_file"]
    score_fs=params["score_fs"]

    scores_dict = dict()
    # Indicate the native structure.
    if decoy_file == ref_file:
        scores_dict["native"] = 1
    else:
        scores_dict["native"] = 0

    scores_dict["pdb_file"] = str(decoy_file)

    # Read in the models.
    m_decoy = IMP.Model()
    h_decoys = IMP.atom.read_multimodel_pdb(
        str(decoy_file),
        m_decoy
    )

    n_states = len(h_decoys)

    m = h_decoys[0].get_model()
    w_p = IMP.Particle(m, "weights")
    w_pid = IMP.isd.Weight.setup_particle(w_p, IMP.algebra.VectorKD([1]*n_states))
    w = IMP.isd.Weight(m, w_pid)

    if not occs:
        occs = list()
        for h in h_decoys:
            occs.append(IMP.atom.Atom(h.get_model(), IMP.atom.Selection(h).get_selected_particle_indexes()[0]).get_occupancy())
    w.set_weights(occs)

    for score_f in score_fs:
        if score_f in ["ff", "bnd", "ang", "dih", "imp", "eps", "nbd"]:
            score = charmm.get_ff_score(
                hs=h_decoys,
                term=score_f
            )
        elif score_f in ["ml", "ls"]:
            # Set f_obs.
            f_obs = miller_ops.get_miller_array(
                f_obs_file=cif_file,
                label="_refln.F_meas_au"
            )
            f_obs = miller_ops.clean_miller_array(f_obs)

            # Set flags.
            status_array = miller_ops.get_miller_array(
                f_obs_file=cif_file,
                label="_refln.status"
            )
            flags = status_array.customized_copy(data=status_array.data()=="f")
            f_obs, flags = f_obs.common_sets(other=flags)

            pids = list()
            for h in h_decoys:
                pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

            results_dict = cctbx_score.get_score(
                hs=h_decoys,
                w=w,
                pids=pids,
                f_obs=f_obs,
                r_free_flags=flags,
                target=score_f,
                u_aniso_file=adp_file
            )
            score = results_dict["score"]
            scores_dict["r_free"] = results_dict["r_free"]
            scores_dict["r_work"] = results_dict["r_work"]
            scores_dict["r_all"] = results_dict["r_all"]
        elif score_f in ["rmsd_avg", "rmsd_ord", "rmsd_dom", "avg_delta_w"]:
            m_ref = IMP.Model()
            h_refs = IMP.atom.read_multimodel_pdb(str(ref_file), m_ref, IMP.atom.AllPDBSelector())

            if score_f == "rmsd_avg":
                f = align_imp.compute_rmsd_between_average
            elif score_f == "rmsd_ord":
                f = align_imp.compute_rmsd_ordered
            elif score_f == "rmsd_dom":
                f = align_imp.compute_rmsd_dom_state
            elif score_f == "avg_delta_w":
                f = align_imp.compute_avg_delta_weight
            score = f(
                    h_0s=h_decoys,
                    h_1s=h_refs,
                    ca_only=True
                )
        elif score_f == "dom_weight":
            hs_ordered = align_imp.get_ordered_hs(h_0s=h_decoys)
            pid_tmp = IMP.atom.Selection(hs_ordered[0]).get_selected_particle_indexes()[0]
            score = IMP.atom.Atom(m, pid_tmp).get_occupancy()
        elif score_f == "w":
            score = list(w.get_weights())
        else:
            score = np.nan

        scores_dict[score_f] = score


    return scores_dict


def score_vs_rmsd(
        params_file,
        pdb_files,
        native_pdb_file,
        native_cif_file,
        flags_file,
        min_res,
        score_fs,
        scores_file
):
    # param_dict = locals()
    # params.write_params(
    #     param_dict=param_dict,
    #     param_file=params_file
    # )

    pool_params = list()

    # Get the number of states.
    for pdb_file in pdb_files:
        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["occs"] = None
        param_dict["ref_file"] = native_pdb_file
        param_dict["cif_file"] = native_cif_file
        param_dict["flags_file"] = flags_file
        param_dict["res"] = min_res
        param_dict["score_fs"] = score_fs
        param_dict["adp_file"] = None

        pool_params.append(param_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_score,
        pool_params
    )

    columns = ["pdb_file", "native"]
    columns.extend(score_fs)
    all_scores_df = pd.DataFrame(columns=columns)
    i = 0
    for score_dict in pool_results:
        print(score_dict)
        for key in score_dict.keys():
            all_scores_df.loc[i, key] = score_dict[key]

        i = i+1

    all_scores_df.to_csv(scores_file)