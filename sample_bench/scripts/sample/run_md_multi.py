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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--cif_files")
    parser.add_argument("--no_k1", action="store_true")
    parser.add_argument("--no_scale", action="store_true")
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--dyn_w_xray", action="store_true")
    parser.add_argument("--start_pdb_file")
    parser.add_argument("--u_aniso_file")
    parser.add_argument("--n_state", type=int)
    parser.add_argument("--init_weights")
    parser.add_argument("--n_cond", type=int)
    parser.add_argument("--ref_pdb_file")
    parser.add_argument("--ref_id", required=False, type=int)
    parser.add_argument("--ref_occs")
    parser.add_argument("--sa")
    parser.add_argument("--steps", type=int)
    parser.add_argument("--bfactor", type=int)
    parser.add_argument("--d_min", type=float)
    args = parser.parse_args()

    params.write_params(vars(args), Path(args.out_dir, "params.txt"))

    ### REPRESENTATION
    pdb_file = Path(args.start_pdb_file)
    n_state = args.n_state
    n_cond = args.n_cond
    pdb_file_n_state = utility.get_n_state_from_pdb_file(pdb_file)

    m, m_0 = IMP.Model(), IMP.Model()
    if pdb_file_n_state == 1:
        hs, h_0s = list(), list()
        for i in range(n_state):
            hs.append(IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector()))
            h_0s.append(IMP.atom.read_pdb(str(pdb_file), m_0, IMP.atom.AllPDBSelector()))
    else:
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        h_0s = IMP.atom.read_multimodel_pdb(str(pdb_file), m_0, IMP.atom.AllPDBSelector())

    # Setup ref_occs.
    # n_ref_state does not have to equal n_state but ref_n_cond has to match n_cond.
    if Path(args.ref_pdb_file).suffix == ".csv":
        ref_df = pd.read_csv(args.ref_pdb_file, index_col=0)
        ref_pdb_file = Path(ref_df.loc[args.ref_id, "pdb"])
    else:
        ref_pdb_file = Path(args.ref_pdb_file)

    ref_n_state = utility.get_n_state_from_pdb_file(ref_pdb_file)
    ref_occs = np.ndarray([ref_n_state, n_cond])

    for cond in range(n_cond):
        for state in range(ref_n_state):
            if Path(args.ref_pdb_file).suffix == ".csv":
                occ = ref_df.loc[args.ref_id, "w_{}_{}".format(state, cond)]
            else:
                occ = args.ref_occs.split(";")[cond].split(",")[state]

            ref_occs[state, cond] = occ
    print("ref_occs: ", ref_occs)

    ref_m = IMP.Model()
    ref_hs = IMP.atom.read_multimodel_pdb(str(ref_pdb_file), ref_m, IMP.atom.AllPDBSelector())

    # Setup occs.
    occs = np.ndarray([n_state, n_cond])
    for cond in range(n_cond):
        if args.init_weights == "rand":
            occs[:, cond] = weights.get_weights(floor=0.05, n_state=n_state)
        elif args.init_weights == "uni":
            occs[:, cond] = [1/n_state]*n_state
        else:
            occs[:, cond] = args.ref_occs.split(";")[cond].split(",")

    print("occs", occs)

    # Setup the weights, indexed by condition.
    ref_ws, ws = list(), list()
    for ws_tmp, occs_tmp, n_state_tmp, m_tmp in [(ref_ws, ref_occs, ref_n_state, ref_m), (ws, occs, n_state, m)]:
        for cond in range(n_cond):
            w_p = IMP.Particle(m_tmp, "weights")
            w_pid = IMP.isd.Weight.setup_particle(w_p, IMP.algebra.VectorKD([1]*n_state_tmp))
            w = IMP.isd.Weight(m_tmp, w_pid)
            w.set_weights(occs_tmp[:,cond])
            ws_tmp.append(w)

    for w in ws:
        w.set_weights_are_optimized(True)

    print("ref_ws:")
    for w in ref_ws:
        print(w.get_weights())

    print("ws:")
    for w in ws:
        print(w.get_weights())

    # Get all particle ids.
    pids = list()
    pids_main_chain = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())
        pids_main_chain.extend(IMP.atom.Selection(h, atom_types=[IMP.atom.AT_CA, IMP.atom.AT_C, IMP.atom.AT_O, IMP.atom.AT_N]).get_selected_particle_indexes())
    pids_side_chain = list(set(pids) - set(pids_main_chain))

    # Setup bfactors.
    if args.bfactor:
        for h in hs:
            for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
                IMP.atom.Atom(m, pid).set_temperature_factor(args.bfactor)

    # Setup waters (if any).
    for h in hs:
        h20s = IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("HET: O  ")).get_selected_particle_indexes()
        for pid in h20s:
            IMP.atom.CHARMMAtom.setup_particle(m, pid, "O")

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
                        if val_str == "A":
                            pids_work = pids
                        elif val_str == "S":
                            pids_work = pids_side_chain
                        elif val_str.isnumeric():
                            pids_work = list()
                            for h in hs:
                                pids_work.extend(IMP.atom.Selection(h, residue_index=int(val_str)).get_selected_particle_indexes())
                        else:
                            raise RuntimeError()

                        val = pids_work

                    sa_step[key] = val

        sa_sched.append(sa_step)

    use_weights = False
    for sa_step in sa_sched:
        if sa_step["w"]:
            use_weights = True

    ps = [m.get_particle(pid) for pid in pids]

    ### SCORING
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
    if args.cif_files:
        cif_files = args.cif_files.split(",")
        pids_xray = pids

        # cif file here is a string.
        for i in range(len(cif_files)):
            cif_file = cif_files[i]
            w = ws[i]

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

            scale = False if args.no_scale else True
            k1 = False if args.no_k1 else True

            r_xray = xray_restraint.XtalRestraint(
                hs=hs,
                pids=pids_xray,
                w=w,
                f_obs=f_obs_array,
                free_flags=flags_array,
                w_xray=args.w_xray,
                update_scale=scale,
                update_k1=k1,
                u_aniso_file=args.u_aniso_file
            )

            rs_xray = IMP.RestraintSet(m, 1.0)
            rs_xray.add_restraint(r_xray)

            r_xrays.append(r_xray)
            rs.append(rs_xray)

            # Setup the weight optimizer state.
            if use_weights and n_state > 1:
                w_os = update_weights_optimizer_state.UpdateWeightsOptimizerState(
                    m=m,
                    hs=hs,
                    w=w,
                    r_xray=r_xray,
                    n_proposals=10,
                    radius=.05
                )

                w_os.set_period(100)
                o_states.append(w_os)

    # Setup the center of mass optimizer state based on the first state of the starting structure.
    pids_ca_0 = list(IMP.atom.Selection(h_0s, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    com_0 = IMP.atom.CenterOfMass.setup_particle(
        IMP.Particle(m_0),
        pids_ca_0
    )

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
        rmsd_all_tracker = trackers.RMSDTracker(
            name="rmsd_{}".format(i),
            hs=hs,
            w=ws[i],
            ref_hs=ref_hs,
            ref_w=ref_ws[i],
            rmsd_func=align_imp.compute_rmsd_between_average,
            ca_only=True
        )
        all_trackers.append(rmsd_all_tracker)

    for i in range(n_cond):
        weight_tracker = trackers.WeightTracker(
            name="w_{}".format(i),
            m=m,
            w=ws[i]
        )
        weight_tracker.set_labels(["w_{}_{}".format(j, i) for j in range(n_state)])
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
        copy_tracker.set_period(10)
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

    for h in hs:
        pids_ca = list(IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
        com_os = com_optimizer_state.CenterOfMassOptimizerState(
            m=m,
            pids=pids_ca,
            com_0=com_0
        )
        o_states.append(com_os)

    if args.steps:
        steps = args.steps
    else:
        steps = 1e99

    molecular_dynamics.molecular_dynamics(
        pdb_dir=tmp_pdb_dir,
        hs=hs,
        r_sets=rs,
        t_step=2,
        n_step=steps,
        sa_sched=sa_sched,
        o_states=o_states
    )

    pdb_tracker.do_evaluate()

    if args.tmp_out_dir:
        copy_tracker.do_evaluate()