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


def pool_score(
        params
):
    scores_dict = score(
        decoy_file=params["decoy_file"],
        occs=params["occs"],
        ref_file=params["ref_file"],
        cif_file=params["cif_file"],
        flags_file=params["flags_file"],
        res=params["res"],
        score_fs=params["score_fs"]
    )

    return scores_dict


def score(
        decoy_file,
        occs,
        ref_file,
        cif_file,
        flags_file,
        res,
        score_fs
):
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

    scores_dict = dict()
    scores_dict["pdb_file"] = str(decoy_file)
    if decoy_file == ref_file:
        scores_dict["native"] = 1
    else:
        scores_dict["native"] = 0

    # Set occupancies.
    if occs:
        if len(occs) != len(h_decoys):
            raise RuntimeError("Length of weights set ({}) is not equal to length of decoy structure set ({}) for decoy ({})".format(len(occs), len(h_decoys), str(decoy_file)))

        for i in range(n_states):
            pids = IMP.atom.Selection(h_decoys[i]).get_selected_particle_indexes()
            for pid in pids:
                IMP.atom.Atom(m_decoy, pid).set_occupancy(occs[i])

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
                hs=h_decoys,
                term=score_term
            )
            scores_dict[score_term] = score
        elif score_term in ["ml", "ls"]:
            # Load the f_obs object.
            f_obs = miller_ops.get_miller_array(
                f_obs_file=cif_file,
                label="_refln.F_meas_au"
            )
            f_obs = miller_ops.clean_miller_array(f_obs)

            f_obs = miller_ops.filter_f_obs_resolution(
                f_obs=f_obs,
                d_max=0,
                d_min=res
            )

            if flags_file:
                status_array = miller_ops.get_miller_array(
                    f_obs_file=flags_file,
                    label="_refln.status"
                )
                flags_array = status_array.customized_copy(data=status_array.data()=="f")
                f_obs, r_free_flags = f_obs.common_sets(other=flags_array)
            else:
                r_free_flags = None

            results_dict = cctbx_score.get_score(
                m=m_decoy,
                f_obs=f_obs,
                r_free_flags=r_free_flags,
                target=score_term
            )
            scores_dict[score_term] = results_dict["score"]
            scores_dict["r_free"] = results_dict["r_free"]
            scores_dict["r_work"] = results_dict["r_work"]
            scores_dict["r_all"] = results_dict["r_all"]
        elif score_term == "rmsd":
            score = align_imp.compute_rmsd_between_average(
                hs_1=h_decoys,
                hs_2=h_refs
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
    if params_file:
        param_dict = locals()
        param_f = open(params_file, "a")
        for key in param_dict.keys():
            print("{:<15}{}\n".format(key, param_dict[key]))
            param_f.write("{:<15}{}\n".format(key, param_dict[key]))
        param_f.write("\n\n")
        param_f.close()

    pool_params = list()

    # Get the number of states.
    m_decoy = IMP.Model()
    h_decoys = IMP.atom.read_multimodel_pdb(
        str(pdb_files[0]),
        m_decoy
    )
    n_states = len(h_decoys)

    for pdb_file in pdb_files:
        param_dict = dict()
        param_dict["decoy_file"] = pdb_file

        if pdb_file == native_pdb_file:
            param_dict["occs"] = [1]
        else:
            param_dict["occs"] = [1/n_states]*n_states

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