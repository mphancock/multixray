from pathlib import Path
import sys
import argparse
import numpy as np
import pandas as pd

import IMP
import IMP.atom
import IMP.core
import IMP.isd
import IMP.algebra

from cctbx.array_family import flex

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
import params
import weight_restraint
import miller_ops
import weights
import utility
import multi_state_multi_condition_model


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--cif_files")
    parser.add_argument("--input_csv")
    parser.add_argument("--job_id", type=int)
    parser.add_argument("--no_k1", action="store_true")
    parser.add_argument("--no_scale", action="store_true")
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--start_pdb_file")
    parser.add_argument("--u_aniso_file")
    parser.add_argument("--n_state", type=int)
    parser.add_argument("--init_weights")
    parser.add_argument("--n_cond", type=int)
    parser.add_argument("--ref_pdb_files")
    parser.add_argument("--ref_id", required=False, type=int)
    parser.add_argument("--ref_occs")
    parser.add_argument("--sa")
    parser.add_argument("--steps", type=int)
    parser.add_argument("--d_min", type=float)
    args = parser.parse_args()

    params.write_params(vars(args), Path(args.out_dir, "params.txt"))

    ## CIF_FILES
    if args.input_csv:
        cif_df = pd.read_csv(Path(args.input_csv), index_col=0)
        cif_files = [Path(cif_file) for cif_file in cif_df.loc[args.job_id, "cifs"].split(",")]
        ref_pdb_files = [Path(ref_file) for ref_file in cif_df.loc[args.job_id, "refs"].split(",")]

        n_cond = len(cif_files)
        ref_n_state = utility.get_n_state_from_pdb_file(ref_pdb_files[0])
        ref_w_mat = np.ndarray(shape=[ref_n_state, n_cond])

        for cond in range(n_cond):
            for state in range(ref_n_state):
                ref_w_mat[state, cond] = cif_df.loc[args.job_id, "ref_occs"].split(";")[cond].split(",")[state]

    # Setup ref models.
    pdb_sel = IMP.atom.NonAlternativePDBSelector()
    ref_msmc_ms = list()
    for cond in range(n_cond):
        ref_msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
            pdb_file=ref_pdb_files[cond],
            w_mat=ref_w_mat
        )

        ref_msmc_ms.append(ref_msmc_m)

    ### REPRESENTATION
    pdb_file = Path(args.start_pdb_file)
    n_state = args.n_state
    pdb_file_n_state = utility.get_n_state_from_pdb_file(pdb_file)

    # # Setup ref_occs.
    # # n_ref_state does not have to equal n_state but ref_n_cond has to match n_cond.
    # if Path(args.ref_pdb_files).suffix == ".csv":
    #     ref_df = pd.read_csv(args.ref_pdb_files, index_col=0)
    #     ref_pdb_files = [Path(ref_df.loc[args.ref_id, "pdb"])]*n_state
    # else:
    #     ref_pdb_files = [Path(ref_pdb_file) for ref_pdb_file in args.ref_pdb_files.split(",")]

    # ref_n_state = utility.get_n_state_from_pdb_file(ref_pdb_files[0])
    # ref_occs = np.ndarray([ref_n_state, n_cond])

    # Setup ref models.
    # ref_msmc_ms = list()
    # for cond in range(n_cond):
    #     ref_m = IMP.Model()
    #     ref_hs = IMP.atom.read_multimodel_pdb(str(ref_pdb_files[cond]), ref_m, pdb_sel)

    #     ref_occs = list()
    #     for state in range(ref_n_state):
    #         if Path(args.ref_pdb_files).suffix == ".csv":
    #             occ = ref_df.loc[args.ref_id, "w_{}_{}".format(state, cond)]
    #         else:
    #             occ = args.ref_occs.split(";")[cond].split(",")[state]

    #         ref_occs.append(occ)

    #     ref_w_mat = np.ndarray(shape=[ref_n_state, n_cond])
    #     ref_w_mat[:, 0] = ref_occs

    #     ref_msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
    #         m=ref_m,
    #         hs=ref_hs,
    #         w_mat=ref_w_mat
    #     )

    #     ref_msmc_ms.append(ref_msmc_m)

    # Setup the multi state multi condition model
    occs = np.ndarray([n_state, n_cond])
    for cond in range(n_cond):
        if args.init_weights == "rand":
            occs[:, cond] = weights.get_weights(floor=0.05, n_state=n_state)
        elif args.init_weights == "uni":
            occs[:, cond] = [1/n_state]*n_state
        else:
            occs[:, cond] = args.ref_occs.split(";")[cond].split(",")

    print("occs", occs)

    msmc_m = multi_state_multi_condition_model.MultiStateMultiConditionModel(
        pdb_file=pdb_file,
        w_mat=occs
    )

    ## SAMPLE
    # Setup simulated annealing schedule.
    sa_sched_str = args.sa.split(";")
    sa_sched_str = [sa_step_str[1:-1] for sa_step_str in sa_sched_str]

    sa_sched = list()
    keys = ["step", "T", "dof", "pdb", "w", "res"]
    for sa_step_str in sa_sched_str:
        sa_step = dict()
        for key_val in sa_step_str.split(","):
            for key in keys:
                if key in key_val:
                    val_str = key_val[len(key):]
                    if key in ["step", "T", "pdb", "w"]:
                        val = int(val_str)
                    elif key == "res":
                        val = float(val_str)
                    else:
                        val = val_str

                    sa_step[key] = val

        sa_sched.append(sa_step)

    use_weights = False
    for sa_step in sa_sched:
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
    rset_charmm.set_weight(1/n_state)
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

            rand_flags = False
            if rand_flags:
                flags_array = f_obs_array.generate_r_free_flags(
                    fraction=0.05,
                    max_free=len(f_obs_array.data())
                )
            else:
                # Set flags from file.
                status_array = miller_ops.get_miller_array(
                    f_obs_file=cif_file,
                    label="_refln.status"
                )
                flags_array = status_array.customized_copy(data=status_array.data()=="f")
            f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

            print("N_OBS: ", len(f_obs_array.data()), len(flags_array.data()))

            com = ref_msmc_ms[i].get_com()
            print("COM: ", com.get_coordinates())

            scale = False if args.no_scale else True
            k1 = False if args.no_k1 else True

            r_xray = xray_restraint.XtalRestraint(
                msmc_m=msmc_m,
                cond=i,
                f_obs=f_obs_array,
                free_flags=flags_array,
                w_xray=args.w_xray,
                # w_xray=args.w_xray/len(cif_files),
                update_scale=scale,
                update_k1=k1,
                u_aniso_file=args.u_aniso_file,
                ref_com=com
            )

            rs_xray = IMP.RestraintSet(m, 1.0)
            rs_xray.add_restraint(r_xray)

            r_xrays.append(r_xray)
            rs.append(rs_xray)

        # Setup the weight optimizer state.
        if use_weights and n_state > 1:
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
    for i in range(len(r_xrays)):
        xray_tracker = trackers.fTracker(
            name="xray_{}".format(i),
            r=r_xrays[i]
        )
        xray_tracker.set_xray_only(True)
        all_trackers.append(xray_tracker)

        r_factor_tracker = trackers.RFactorTracker(
            name="r_factor_{}".format(i),
            r_xray=r_xrays[i]
        )
        r_factor_tracker.set_labels(["r_free_{}".format(i), "r_work_{}".format(i)])
        all_trackers.append(r_factor_tracker)

    for i in range(n_cond):
        ref_msmc_m = ref_msmc_ms[i]
        rmsd_all_tracker = trackers.RMSDTracker(
            name="rmsd_{}".format(i),
            rmsd_func=align_imp.compute_rmsd_between_average,
            hs_0=msmc_m.get_hs(),
            hs_1=ref_msmc_m.get_hs(),
            pids_0=msmc_m.get_ca_pids(0),
            pids_1=ref_msmc_m.get_ca_pids(0),
            occs_0=msmc_m.get_w_mat()[:, i],
            occs_1=ref_msmc_m.get_w_mat()[:, 0]
        )
        all_trackers.append(rmsd_all_tracker)

    weight_tracker = trackers.WeightMatTracker(
        name="w",
        msmc_m=msmc_m
    )
    all_trackers.append(weight_tracker)

    # pdb writer trackers.
    pdb_dir = Path(args.out_dir, "pdbs")
    pdb_dir.mkdir()

    if args.tmp_out_dir:
        tmp_pdb_dir = Path(args.tmp_out_dir, "pdbs")
        tmp_pdb_dir.mkdir()
    else:
        tmp_pdb_dir = pdb_dir

    pdb_tracker = pdb_writer.PDBWriterTracker(
        name="pdb",
        hs=hs,
        pdb_dir=tmp_pdb_dir,
        log_pdb_dir=pdb_dir
    )
    pdb_tracker.set_period(10)
    all_trackers.append(pdb_tracker)

    # This freq needs to be equal to the logging frequency to ensure there are never more log entries than corresponding pdb files.
    if args.tmp_out_dir:
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
        log_file=Path(args.out_dir, "log.csv"),
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

    if args.steps:
        steps = args.steps
    else:
        steps = 1e99

    molecular_dynamics.molecular_dynamics(
        msmc_m=msmc_m,
        r_sets=rs,
        t_step=2,
        n_step=steps,
        sa_sched=sa_sched,
        o_states=o_states
    )

    pdb_tracker.do_evaluate()

    if args.tmp_out_dir:
        copy_tracker.do_evaluate()