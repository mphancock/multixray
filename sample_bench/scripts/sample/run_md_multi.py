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

from cctbx.crystal import symmetry

sys.path.append(str(Path(Path.home(), "xray/src")))
from charmm import get_charmm_restraint_set, CHARMMDerivHolder
import xray_restraint
import trackers
import log_statistics
import align_imp
import pdb_writer
import com_optimizer_state
import update_weights_optimizer_state
from simulated_annealing import SimulatedAnnealing, SimulatedAnnealingSchedule
from params import write_params_txt, write_params_csv, read_job_csv
from miller_ops import get_crystal_symmetry, get_f_obs, get_flags, clean_miller_array
import weights
import utility
import multi_state_multi_condition_model
from derivatives import evaluate_df_dict, DerivScoreState


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--job_csv_file")
    parser.add_argument("--job_id")
    parser.add_argument("--write", action="store_true")
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
    w_xray = params_dict["w_xray"]

    # write_params_csv(param_dict=params_dict, param_file=Path(out_dir, "params.csv"))

    # Setup the multi state multi condition model
    crystal_symmetries = list()

    for cond in range(J):
        ## may need to create a fake crystal symmetry if there are no cif files
        if not cif_files[cond]:
            crystal_symmetry = symmetry([100, 100, 100, 90, 90, 90], "P1")
            crystal_symmetries.append(crystal_symmetry)
        else:
            crystal_symmetries.append(get_crystal_symmetry(
                cif_file=cif_files[cond]
            ))

    ## REPRESENTATION
    msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
        pdb_files=params_dict["init_pdbs"],
        w_mat=params_dict["init_w_mat"],
        crystal_symmetries=crystal_symmetries
    )

    ref_msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
        pdb_files=params_dict["ref_pdbs"],
        w_mat=params_dict["ref_w_mat"],
        crystal_symmetries=crystal_symmetries
    )

    ### SCORING
    m, hs = msmc_m.get_m(), msmc_m.get_hs()
    r_charmm = get_charmm_restraint_set(m, hs)

    if params_dict["refine"]:
        print("REFINING STARTING MODEL")

        sf_refine = IMP.core.RestraintsScoringFunction([r_charmm])
        cg = IMP.core.ConjugateGradients(msmc_m.get_m())
        cg.set_scoring_function(sf_refine)
        cg.optimize(1000)

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

        # com = ref_msmc_ms[i].get_com()

        xray_freq = params_dict["xray_freq"]
        r_xray = xray_restraint.XtalRestraint(
            msmc_m=msmc_m,
            cond=i,
            f_obs=f_obs,
            free_flags=flags,
            w_xray=w_xray,
            update_scale=True,
            update_k1=True,
            remove_outliers=True,
            update_freq=xray_freq,
            charmm_holder=charmm_deriv_holder,
            ref_com=None
        )

        rset_xray.add_restraint(r_xray)

    r_xrays = [rset_xray.get_restraint(i) for i in range(rset_xray.get_number_of_restraints())]

    # Setup the weight optimizer state.
    w_os = update_weights_optimizer_state.UpdateWeightsOptimizerState(
        msmc_m=msmc_m,
        r_xrays=r_xrays,
        n_proposals=10,
        radius=.05,
        write=args.write
    )
    w_os.set_period(100)

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
        r=r_charmm
    )
    all_trackers.append(ff_tracker)

    # Add the trackers for each xtal restraint.
    ## if there is 1 condition there may be 1 or 0 xray restraints
    ## else there must be an extra restraint for each condition
    n_xray_rs = rset_xray.get_number_of_restraints()
    assert n_xray_rs == 0 or n_xray_rs == J

    for cond in range(n_xray_rs):
        r_xray = rset_xray.get_restraint(cond)
        if not cif_files[cond]:
            continue

        cif_name = cif_files[cond].stem
        xray_tracker = trackers.fTracker(
            name="xray_{}".format(cif_name),
            r=r_xray
        )
        xray_tracker.set_xray_only(True)
        all_trackers.append(xray_tracker)

        r_factor_tracker = trackers.RFactorTracker(
            name="r_factor_{}".format(cif_name),
            r_xray=r_xray,
            labels=["r_free_{}".format(cif_name), "r_work_{}".format(cif_name)]
        )
        all_trackers.append(r_factor_tracker)

    rmsd_tracker = trackers.MultiStateMultiConditionRMSDTracker(
        name="rmsd",
        msmc_0=msmc_m,
        msmc_1=ref_msmc_m
    )
    all_trackers.append(rmsd_tracker)

    weight_labels = list()
    for state in range(N):
        for cond in range(J):
            if cif_files[cond]:
                weight_name = "w_{}_{}".format(state, cif_files[cond].stem)
            else:
                weight_name = "w_{}_{}".format(state, cond)

            weight_labels.append(weight_name)

    weight_tracker = trackers.WeightMatTracker(
        name="w",
        msmc_m=msmc_m,
        labels=weight_labels
    )
    all_trackers.append(weight_tracker)

    pdb_tracker = pdb_writer.PDBWriterTracker(
        name="pdb",
        msmc_m=msmc_m,
        pdb_dir=tmp_pdb_dir,
        log_pdb_dir=pdb_dir
    )
    pdb_tracker.set_period(10)
    all_trackers.append(pdb_tracker)

    # This freq needs to be equal to the logging frequency to ensure there are never more log entries than corresponding pdb files.
    if tmp_out_dir:
        copy_tracker = trackers.CopyTracker(
            name="copy",
            m=m,
            source_dir=tmp_out_dir,
            dest_dir=out_dir
        )
        copy_tracker.set_xray_only(False)
        copy_tracker.set_period(100)
        all_trackers.append(copy_tracker)

    ## magnitude trackers
    ## average velocity magnitude
    velocity_mag_tracker = trackers.VelocityMagnitudeTracker(
        name="vel_mag",
        m=m,
        pids=msmc_m.get_pids()
    )
    all_trackers.append(velocity_mag_tracker)

    ## average charmmm derivative magnitude
    dcharmm_mag_tracker = trackers.dfMagnitudeTracker(
        name="dcharmm_mag",
        m=m,
        pids=msmc_m.get_pids(),
        r=r_charmm,
        scale=1
    )
    all_trackers.append(dcharmm_mag_tracker)

    ## average xray derivative magnitude
    if len(r_xrays) > 0:
        for i in range(rset_xray.get_number_of_restraints()):
            cif_name = cif_files[i].stem
            dxray_magnitude_tracker = trackers.dfMagnitudeTracker(
                name="dxray_{}_mag".format(cif_name),
                m=m,
                pids=msmc_m.get_pids(),
                r=rset_xray.get_restraint(i),
                scale=1000000
            )
            all_trackers.append(dxray_magnitude_tracker)

            ## add one for each state
            for state in range(N):
                dxray_state_mag_tracker = trackers.dfMagnitudeTracker(
                    name="dxray_{}_mag_{}".format(cif_name, state),
                    m=m,
                    pids=msmc_m.get_ca_pids_in_state(state),
                    r=rset_xray.get_restraint(i),
                    scale=1000000
                )
                all_trackers.append(dxray_state_mag_tracker)

    ## average derivative magnitude
    dXYZ_mag_tracker = trackers.dXYZMagnitudeTracker(
        name="dxyz_mag",
        m=m,
        pids=msmc_m.get_pids()
    )

    ## track the weight of the xray restraints
    if len(r_xrays) > 0:
        wxray_tracker = trackers.XrayWeightTracker(
            name="wxray",
            m=msmc_m.get_m(),
            r_xray=rset_xray.get_restraint(0)
        )
        all_trackers.append(wxray_tracker)

    # track the temperature but turn off until md is added
    temp_tracker = trackers.TempTracker(
        name="temp",
        m=m
    )
    temp_tracker.set_off()
    all_trackers.append(temp_tracker)

    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=Path(tmp_out_dir, "log.csv"),
        log_freq=1,
        write=args.write
    )

    r_charmm.evaluate(calc_derivs=True)
    rset_xray.evaluate(calc_derivs=True)

    log_ostate.update()

    ## setup score state for updating derivatives
    deriv_score_state = DerivScoreState(
        m=m,
        pids=msmc_m.get_pids(),
        charmm_holder=charmm_deriv_holder,
        xray_rs=r_xrays,
        w_xray=w_xray
    )
    m.add_score_state(deriv_score_state)

    # Need one absolute center of mass
    # com_os = com_optimizer_state.CenterOfMassOptimizerState(
    #     m=m,
    #     pids=msmc_m.get_ca_pids(),
    #     ref_com=ref_msmc_ms[0].get_com()
    # )
    com_os = None

    ## Max step of 2 just to run 1 cycle of the sampling schedule
    sa_sched = SimulatedAnnealingSchedule(sa_string=params_dict["sample_sched_str"])
    sa = SimulatedAnnealing(
        msmc_m=msmc_m,
        rset_xray=rset_xray,
        r_charmm=r_charmm,
        n_step=2,
        sa_sched=sa_sched,
        log_o_state=log_ostate,
        weight_o_state=w_os,
        com_o_state=com_os,
        weight_thermo=params_dict["weight_thermo"],
        vel_thermo=params_dict["vel_thermo"]
    )
    sa.run()

    pdb_tracker.do_evaluate()

    if tmp_out_dir:
        copy_tracker.do_evaluate()