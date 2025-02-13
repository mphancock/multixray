import sys
from pathlib import Path
import shutil
import numpy as np
import math

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import params
import trackers
import log_statistics
import pdb_writer
from utility import pool_read_pdb
from derivatives import DerivScoreState


def read_pdb_and_refine_to_max_ff(
    pool_params
):
    pdb_file = pool_params["pdb_file"]
    out_pdb_file = pool_params["out_pdb_file"]
    max_ff = pool_params["max_ff"]
    log_file = pool_params["log_file"]

    if pdb_file.exists():
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    else:
        return RuntimeError("{} file does not exist.".format(pdb_file))

    refine_hs_to_max_ff(hs=hs, max_ff=max_ff, log_file=log_file)

    if out_pdb_file:
        IMP.atom.write_multimodel_pdb(hs, str(out_pdb_file))

    return out_pdb_file


def refine_hs_to_max_ff(
    hs,
    max_ff,
    log_file
):
    for i in range(len(hs)):
        print(i)
        h = hs[i]
        m = h.get_model()

        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        # print(len(pids))
        for pid in pids:
            IMP.core.XYZ(m, pid).set_coordinates_are_optimized(True)

        ps = [m.get_particle(pid) for pid in pids]

        rs = list()
        rset_charmm = IMP.RestraintSet(m, 1.0)
        charmm_rs = charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
        rset_charmm.add_restraints(charmm_rs)
        rset_charmm.set_weight(1)
        rs.append(rset_charmm)

        # If the initial structure is non-physiological, don't refine
        ff_cur = rset_charmm.evaluate(False)
        if math.isnan(ff_cur) or ff_cur > 5000000:
            keep_refining = False
        else:
            keep_refining = True

        while keep_refining:
            # ff_cur = rset_charmm.evaluate(False)
            # print(ff_cur)
            refine_h(
                h=h,
                rs=rs,
                n_step=10,
                log_file=log_file
            )

            ff_new = rset_charmm.evaluate(False)
            print(ff_cur, ff_new)

            if ff_new < max_ff:
                keep_refining = False
            elif abs(ff_cur - ff_new) < 100:
                keep_refining = False
            elif ff_new > 2500000:
                keep_refining = False

            ff_cur = ff_new


"""
I read and refine in a single function to reduce memory usage versus loading all hierarchies/models at once.
"""
def read_pdb_and_refine(
    pool_params
):
    pdb_file = pool_params["pdb_file"]
    out_pdb_file = pool_params["out_pdb_file"]
    n_step = pool_params["n_step"]
    max_ff = pool_params["max_ff"]
    log_file = pool_params["log_file"]

    if pdb_file.exists():
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    else:
        return RuntimeError("{} file does not exist.".format(pdb_file))

    if n_step:
        refine_hs(hs=hs, n_step=n_step, log_file=log_file)

    if out_pdb_file:
        IMP.atom.write_multimodel_pdb(hs, str(out_pdb_file))

    return out_pdb_file


def pool_refine(
    pool_param
):
    hs = pool_param["hs"]
    n_step = pool_param["n_step"]
    log_file = pool_param["log_file"]

    return refine_hs(hs=hs, n_step=n_step, log_file=log_file)


def refine_hs(
    hs,
    n_step,
    log_file
):
    for i in range(len(hs)):
        h = hs[i]

        refine_h(
            h=h,
            n_step=n_step,
            log_file=log_file
        )

    return m, hs


def refine_h(
    h,
    rs,
    n_step
):
    m = h.get_model()

    sf = IMP.core.RestraintsScoringFunction(rs)
    cg = IMP.core.ConjugateGradients(m)
    cg.set_scoring_function(sf)

    for o_state in o_states:
        cg.add_optimizer_state(o_state)

    cg.optimize(n_step)


sys.path.append(str(Path(Path.home(), "xray/src")))
from params import write_params_txt, write_params_csv, read_job_csv
from miller_ops import get_crystal_symmetry, get_f_obs, get_flags, clean_miller_array
import multi_state_multi_condition_model
from charmm import get_charmm_restraint_set, CHARMMDerivHolder
import xray_restraint


def refine_posterior(
    pdb_file,
    cif_files,
    w_mat,
    ref_pdb_file,
    n_step
):
    N, J = w_mat.shape

    crystal_symmetries = []
    for cond in range(J):
        crystal_symmetries.append(get_crystal_symmetry(
            cif_file=cif_files[cond]
        ))

    msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
        pdb_files=[pdb_file],
        w_mat=w_mat,
        crystal_symmetries=crystal_symmetries
    )

    ### SCORING
    m, hs = msmc_m.get_m(), msmc_m.get_hs()
    r_charmm = get_charmm_restraint_set(m, hs)

    # cif file here is a string.
    rset_xray = IMP.RestraintSet(m, 1.0)
    charmm_deriv_holder = CHARMMDerivHolder()
    for i in range(len(cif_files)):
        cif_file = cif_files[i]
        if not cif_file:
            continue

        f_obs = get_f_obs(cif_file=cif_file)
        f_obs = clean_miller_array(f_obs)

        flags = get_flags(cif_file=cif_file)
        f_obs, flags = f_obs.common_sets(other=flags)

        r_xray = xray_restraint.XtalRestraint(
            msmc_m=msmc_m,
            cond=i,
            f_obs=f_obs,
            free_flags=flags,
            w_xray=1/len(cif_files),
            update_scale=True,
            update_k1=True,
            remove_outliers=False,
            update_freq=1,
            charmm_holder=charmm_deriv_holder,
            ref_com=None
        )

        rset_xray.add_restraint(r_xray)


    r_xrays = [rset_xray.get_restraint(i) for i in range(rset_xray.get_number_of_restraints())]

    for r_xray in r_xrays:
        r_xray.evaluate(True)

    ## setup score state for updating derivatives
    deriv_score_state = DerivScoreState(
        m=m,
        pids=msmc_m.get_pids(),
        charmm_holder=charmm_deriv_holder,
        xray_rs=r_xrays,
        w_xray=1
    )
    m.add_score_state(deriv_score_state)

    sf_refine = IMP.core.RestraintsScoringFunction([r_charmm, rset_xray])
    cg = IMP.core.ConjugateGradients(msmc_m.get_m())
    cg.set_scoring_function(sf_refine)
    cg.optimize(n_step)

    IMP.atom.write_multimodel_pdb(hs, str(ref_pdb_file))

    score_dict = dict()
    score_dict["ff"] = r_charmm.evaluate(True)

    for i in range(J):
        cif_name = cif_files[i].stem
        score_dict["r_free_{}".format(cif_name)] = r_xrays[i].get_r_free()
        score_dict["r_work_{}".format(cif_name)] = r_xrays[i].get_r_work()
        score_dict["xray_{}".format(cif_name)] = r_xrays[i].get_score()

    return score_dict
