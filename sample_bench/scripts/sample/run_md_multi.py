from pathlib import Path
import sys
import argparse
import numpy as np
import pandas as pd
import random
import math
import shutil

import IMP
import IMP.atom
import IMP.core
import IMP.isd
import IMP.algebra

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import xray_restraint
import trackers
import com_restraint
import log_statistics
import align_imp
import pdb_writer
import com_optimizer_state
import update_weights_optimizer_state
import reset
import molecular_dynamics
from params import write_params_txt, write_params_pickle, write_params_csv, read_job_csv
import weight_restraint
import miller_ops
import weights
import utility
from utility import get_cif_and_ref_from_input_df, get_sa_sched_from_string
import multi_state_multi_condition_model


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--job_csv_file")
    parser.add_argument("--job_id")
    args = parser.parse_args()

    ## Be able to read from a params csv to make rerunning jobs easier
    job_csv_file = Path(args.job_csv_file)
    params_dict = read_job_csv(job_csv_file=job_csv_file, job_id=int(args.job_id))
    print(params_dict)

    ## Always have to pass out_dir and tmp_out_dir to the script
    out_dir = Path(args.out_dir)
    shutil.rmtree(out_dir, ignore_errors=True)
    pdb_dir = Path(out_dir, "pdbs")
    pdb_dir.mkdir(parents=True, exist_ok=True)

    tmp_out_dir = Path(args.tmp_out_dir)
    tmp_pdb_dir = Path(tmp_out_dir, "pdbs")
    tmp_pdb_dir.mkdir(parents=True, exist_ok=True)

    N = params_dict["N"]
    J = params_dict["J"]
    cif_files = params_dict["cifs"]
    ref_pdb_files = params_dict["refs"]
    ref_w_mat = params_dict["ref_w_mat"]
    w_xray = params_dict["w_xray"]
    sample_sched = params_dict["sample_sched"]

    # Optional params
    start_pdb_file = params_dict["start_pdb_file"]
    init_weights = params_dict["init_weights"]

    write_params_csv(param_dict=params_dict, param_file=Path(out_dir, "params.csv"))

    # Setup ref models.
    pdb_sel = IMP.atom.NonAlternativePDBSelector()
    ref_msmc_ms = list()
    for cond in range(J):
        ref_msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
            pdb_file=ref_pdb_files[cond],
            w_mat=ref_w_mat
        )

        ref_msmc_ms.append(ref_msmc_m)

    ### REPRESENTATION
    if start_pdb_file:
        pdb_file = Path(start_pdb_file)
    else:
        pdb_file = random.choice(ref_pdb_files)

    # Setup the multi state multi condition model
    w_mat = np.ndarray([N, J])
    for cond in range(J):
        if init_weights:
            w_mat[:, cond] = ref_occs.split(";")[cond].split(",")
        else:
            w_mat[:, cond] = weights.get_weights(floor=0.05, n_state=N)

    msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
        pdb_file=pdb_file,
        w_mat=w_mat
    )

    ## SAMPLE
    # Setup simulated annealing schedule.

    use_weights = False
    for sa_step in sample_sched:
        if sa_step["w"]:
            use_weights = True

    ### SCORING
    m, hs = msmc_m.get_m(), msmc_m.get_hs()
    rs = list()
    rset_charmm = IMP.RestraintSet(m, 1.0)
    for h in hs:
        charmm_rs = charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
        rset_charmm.add_restraints(charmm_rs)

    # Is turning this off going to be really bad?
    # rset_charmm.set_weight(1/N)
    rs.append(rset_charmm)

    # List of the atom and weight restraints.
    r_xrays, wrs = list(), list()
    o_states = list()
    if cif_files:
        # cif file here is a string.
        for i in range(len(cif_files)):
            cif_file = cif_files[i]

            f_obs_array = miller_ops.get_miller_array(
                f_obs_file=cif_file,
                label="_refln.F_meas_au"
            )
            f_obs_array = miller_ops.clean_miller_array(f_obs_array)

            # Set flags from file.
            status_array = miller_ops.get_miller_array(
                f_obs_file=cif_file,
                label="_refln.status"
            )
            flags_array = status_array.customized_copy(data=status_array.data()=="f")
            f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

            # print("N_OBS: ", len(f_obs_array.data()), len(flags_array.data()))

            com = ref_msmc_ms[i].get_com()

            r_xray = xray_restraint.XtalRestraint(
                msmc_m=msmc_m,
                cond=i,
                f_obs=f_obs_array,
                free_flags=flags_array,
                w_xray=w_xray,
                # w_xray=args.w_xray/len(cif_files),
                update_scale=True,
                update_k1=True,
                u_aniso_file=None,
                ref_com=com
            )

            rs_xray = IMP.RestraintSet(m, 1.0)
            rs_xray.add_restraint(r_xray)

            r_xrays.append(r_xray)
            rs.append(rs_xray)

        # Setup the weight optimizer state.
        if use_weights and N > 1:
            w_os = update_weights_optimizer_state.UpdateWeightsOptimizerState(
                msmc_m=msmc_m,
                r_xrays=r_xrays,
                n_proposals=10,
                radius=.05
            )

            w_os.set_period(100)
            o_states.append(w_os)

    all_trackers = list()
    step_tracker = trackers.StepTracker(
        name="step",
        m=m
    )
    all_trackers.append(step_tracker)

    time_tracker = trackers.TimeTracker(
        name="time",
        m=m
    )
    all_trackers.append(time_tracker)

    ff_tracker = trackers.fTracker(
        name="ff",
        r=rset_charmm
    )
    all_trackers.append(ff_tracker)

    # Add the trackers for each xtal restraint.
    for cond in range(J):
        cif_name = cif_files[cond].stem
        xray_tracker = trackers.fTracker(
            name="xray_{}".format(cif_name),
            r=r_xrays[cond]
        )
        xray_tracker.set_xray_only(True)
        all_trackers.append(xray_tracker)

        r_factor_tracker = trackers.RFactorTracker(
            name="r_factor_{}".format(cif_name),
            r_xray=r_xrays[cond],
            labels=["r_free_{}".format(cif_name), "r_work_{}".format(cif_name)]
        )
        all_trackers.append(r_factor_tracker)

    for cond in range(J):
        cif_name = cif_files[cond].stem
        ref_msmc_m = ref_msmc_ms[cond]
        rmsd_all_tracker = trackers.RMSDTracker(
            name="rmsd_{}".format(cif_name),
            rmsd_func=align_imp.compute_rmsd_between_average,
            hs_0=msmc_m.get_hs(),
            hs_1=ref_msmc_m.get_hs(),
            pids_0=msmc_m.get_ca_pids(0),
            pids_1=ref_msmc_m.get_ca_pids(0),
            occs_0=msmc_m.get_w_mat()[:, cond],
            occs_1=ref_msmc_m.get_w_mat()[:, 0]
        )
        all_trackers.append(rmsd_all_tracker)

    weight_labels = list()
    for state in range(N):
        for cond in range(J):
            weight_labels.append("w_{}_{}".format(state, cif_files[cond].stem))

    weight_tracker = trackers.WeightMatTracker(
        name="w",
        msmc_m=msmc_m,
        labels=weight_labels
    )
    all_trackers.append(weight_tracker)

    pdb_tracker = pdb_writer.PDBWriterTracker(
        name="pdb",
        hs=hs,
        pdb_dir=tmp_pdb_dir,
        log_pdb_dir=pdb_dir
    )
    pdb_tracker.set_period(10)
    all_trackers.append(pdb_tracker)

    # This freq needs to be equal to the logging frequency to ensure there are never more log entries than corresponding pdb files.
    if tmp_out_dir:
        copy_tracker = pdb_writer.PDBCopyTracker(
            name="copy",
            m=m,
            source_dir=tmp_pdb_dir,
            dest_dir=pdb_dir
        )
        copy_tracker.set_xray_only(False)
        copy_tracker.set_period(100)
        all_trackers.append(copy_tracker)

    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=Path(out_dir, "log.csv"),
        log_freq=100
    )

    for r in rs:
        r.evaluate(calc_derivs=True)
    log_ostate.update()

    o_states.append(log_ostate)

    # Need one absolute center of mass
    com_os = com_optimizer_state.CenterOfMassOptimizerState(
        m=m,
        pids=msmc_m.get_all_ca_pids(),
        ref_com=ref_msmc_ms[0].get_com()
    )
    o_states.append(com_os)

    ## Max step of 2 just to run 1 cycle of the sampling schedule
    molecular_dynamics.molecular_dynamics(
        msmc_m=msmc_m,
        r_sets=rs,
        t_step=2,
        n_step=2,
        sa_sched=sample_sched,
        o_states=o_states
    )

    pdb_tracker.do_evaluate()

    if tmp_out_dir:
        copy_tracker.do_evaluate()