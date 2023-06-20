import sys
from pathlib import Path
import multiprocessing
import itertools
import pandas as pd
import numpy as np

import IMP, IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp
# import cctbx_scores
import cctbx_score
import charmm
import miller_ops
import params


def get_params(
        params_dict,
        combos
):
    if combos:
        itertools_args = list()
        for param_name in list(params_dict.keys()):
            itertools_args.append(params_dict[param_name])

        params = list(itertools.product(*tuple(itertools_args)))

    else:
        params = list()
        for i in range(len(params_dict["decoy"])):
            param_entry = list()
            for param_name in list(params_dict.keys()):
                param_entry.append(params_dict[param_name][i])

            params.append(tuple(param_entry))

    return params


"""
The function "pool_score" takes in a dictionary of parameters (params) and calculates scores based on the provided information. It reads in decoy and reference files in PDB format, sets occupancies, retrieves Miller array data, sets flags, and calls the score function to calculate scores. The scores are stored in a dictionary (scores_dict) and returned.

Parameters:
    params: A dictionary containing various parameters required for score calculation, including "decoy_file" (path to decoy file), "occs" (occupancy values), "ref_file" (path to reference file), "cif_file" (path to CIF file), "res" (resolution), and "score_fs" (list of score terms).

Return:
    scores_dict: A dictionary containing the calculated scores.

"""
def pool_score(
        params
):
    decoy_file=params["decoy_file"]
    occs=params["occs"]
    ref_file=params["ref_file"]
    cif_file=params["cif_file"]
    res=params["res"]
    score_fs=params["score_fs"]

    scores_dict = dict()
    # Indicate the native structure.
    if decoy_file == ref_file:
        scores_dict["native"] = 1
    else:
        scores_dict["native"] = 0

    scores_dict["pdb_file"] = str(decoy_file)

    # Read in the models.
    m_ref = IMP.Model()
    h_refs = IMP.atom.read_multimodel_pdb(
        str(ref_file),
        m_ref
    )

    m_decoy = IMP.Model()
    h_decoys = IMP.atom.read_multimodel_pdb(
        str(decoy_file),
        m_decoy
    )

    n_states = len(h_decoys)

    # Set occupancies. -1 means to use uniform occupancies.
    if not occs:
        occs = [1/n_states]*n_states

    if len(occs) != len(h_decoys):
        raise RuntimeError("Length of weights set ({}) is not equal to length of decoy structure set ({}) for decoy ({})".format(len(occs), len(h_decoys), str(decoy_file)))

    for i in range(n_states):
        pids = IMP.atom.Selection(h_decoys[i]).get_selected_particle_indexes()
        for pid in pids:
            IMP.atom.Atom(m_decoy, pid).set_occupancy(occs[i])

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

    score_value_dict = score(
        hs=h_decoys,
        hs_0=h_refs,
        f_obs=f_obs,
        flags=flags,
        res=res,
        score_fs=score_fs
    )
    scores_dict.update(score_value_dict)

    return scores_dict


"""
The function "score" takes in various parameters such as hs, hs_0, cif_file, flags_file, res, and score_fs, and calculates different scores based on the given parameters. It uses the CHARMM program to compute force field (ff) scores and the CCTBX library to calculate maximum likelihood (ml) and least squares (ls) scores. The function also computes the root-mean-square deviation (rmsd) between two sets of input data. The scores are stored in a dictionary with the corresponding score term as the key.

Parameters:
    hs: A list of hierarchies. The occupancies of these heirachies will be used in score calculations.
    hs_0: A list of reference hierachies used for computing RMSD.
    f_obs: A miller array containing the observed structure factor amplitudes.
    flags: A miller array containing the r_free flags for the observed structure factor amplitudes.
    res: A numerical value representing the resolution. *****NOT IMPLEMENTED*****
    score_fs: A list of score terms for which scores need to be calculated.

Return:
    scores_dict: A dictionary containing the calculated scores, with score terms as keys.
"""
def score(
        hs,
        hs_0,
        f_obs,
        flags,
        res,
        score_fs
):
    m = hs[0].get_model()
    scores_dict = dict()

    if "total" in score_fs:
        score_terms = score_fs.copy()
        score_terms.remove("total")

        if "ml" not in score_terms:
            score_terms.append("ml")
        if "ff" not in score_terms:
            score_terms.append("ff")
    else:
        score_terms = score_fs

    for score_term in score_terms:
        # First evaluate all ff terms.
        if score_term in ["ff", "bnd", "ang", "dih", "imp", "eps", "nbd"]:
            score = charmm.get_ff_score(
                hs=hs,
                term=score_term
            )
            scores_dict[score_term] = score
        elif score_term in ["ml", "ls"]:
            results_dict = cctbx_score.get_score(
                m=m,
                f_obs=f_obs,
                r_free_flags=flags,
                target=score_term
            )
            scores_dict[score_term] = results_dict["score"]
            scores_dict["r_free"] = results_dict["r_free"]
            scores_dict["r_work"] = results_dict["r_work"]
            scores_dict["r_all"] = results_dict["r_all"]
        elif score_term == "rmsd":
            score = align_imp.compute_rmsd_between_average(
                hs_1=hs,
                hs_2=hs_0
            )
            scores_dict["rmsd"] = score
        else:
            scores_dict[score_term] = np.nan

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
    param_dict = locals()
    params.write_params(
        param_dict=param_dict,
        param_file=params_file
    )

    pool_params = list()

    # Get the number of states.
    m_decoy = IMP.Model()
    h_decoys = IMP.atom.read_multimodel_pdb(
        str(pdb_files[0]),
        m_decoy
    )

    for pdb_file in pdb_files:
        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["occs"] = None
        param_dict["ref_file"] = native_pdb_file
        param_dict["cif_file"] = native_cif_file
        param_dict["flags_file"] = flags_file
        param_dict["res"] = min_res
        param_dict["score_fs"] = score_fs

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