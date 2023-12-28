from pathlib import Path
import sys
import argparse
import shutil
import random
import numpy as np

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


def get_n_state_from_pdb_file(pdb_file):
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    return len(hs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--cif_files")
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--dyn_w_xray", action="store_true")
    parser.add_argument("--start_pdb_file")
    parser.add_argument("--u_aniso_file")
    parser.add_argument("--n_state")
    parser.add_argument("--init_weights")
    parser.add_argument("--ref_pdb_files")
    parser.add_argument("--sa")
    parser.add_argument("--steps", type=int)
    parser.add_argument("--bfactor", type=int)
    parser.add_argument("--dropout", action="store_true")
    parser.add_argument("--d_min", type=float)
    parser.add_argument("--rand_noise", action="store_true")
    parser.add_argument("--main_chain", action="store_true")
    args = parser.parse_args()

    params.write_params(vars(args), Path(args.out_dir, "params.txt"))

    ### REPRESENTATION
    pdb_file = Path(args.start_pdb_file)
    n_states = int(args.n_state)
    pdb_file_n_state = get_n_state_from_pdb_file(pdb_file)

    m, m_0 = IMP.Model(), IMP.Model()
    if pdb_file_n_state == 1:
        hs, h_0s = list(), list()
        for i in range(n_states):
            hs.append(IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector()))
            h_0s.append(IMP.atom.read_pdb(str(pdb_file), m_0, IMP.atom.AllPDBSelector()))
    else:
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        h_0s = IMP.atom.read_multimodel_pdb(str(pdb_file), m_0, IMP.atom.AllPDBSelector())

    # There are multiple ref pdb files because there should be one for each Xray dataset.
    ref_pdb_files = [Path(ref_pdb_file) for ref_pdb_file in args.ref_pdb_files.split(",")]
    all_ref_ms, all_ref_hs = list(), list()
    for ref_pdb_file in ref_pdb_files:
        ref_m = IMP.Model()
        ref_hs = IMP.atom.read_multimodel_pdb(str(ref_pdb_file), ref_m, IMP.atom.AllPDBSelector())
        all_ref_ms.append(ref_m)
        all_ref_hs.append(ref_hs)

    ws = list()
    if args.cif_files:
        n_cif = len(args.cif_files.split(","))
        if len(all_ref_hs) != n_cif:
            raise RuntimeError("Number of cif files ({}) does not match number of ref pdb files ({})".format(n_cif, len(all_ref_hs)))

        for i in range(n_cif):
            # Setup the weights.
            w_p = IMP.Particle(m, "weights")
            w_pid = IMP.isd.Weight.setup_particle(w_p, IMP.algebra.VectorKD([1]*n_states))
            w = IMP.isd.Weight(m, w_pid)

            if args.init_weights == "ref":
                w_vals = weights.get_weights_from_hs(all_ref_hs[i])
            elif args.init_weights == "rand":
                w_vals = weights.get_weights(
                    floor=0.05,
                    ws_cur=[1/n_states]*n_states,
                    sigma=None
                )
            elif args.init_weights == "uni":
                weights = [1/n_states]*n_states
            else:
                init_weights = args.init_weights.split(",")
                init_weights = [float(w) for w in init_weights]
                if len(init_weights) != n_states:
                    raise RuntimeError("Number of initial weights does not match number of states.")

            w.set_weights(w_vals)
            w.set_weights_are_optimized(True)
            print(w.get_weights())

            ws.append(w)

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
    rset_charmm.set_weight(1/n_states)
    rs.append(rset_charmm)

    # List of the atom and weight restraints.
    r_xrays, wrs = list(), list()
    o_states = list()
    if args.cif_files:
        cif_files = args.cif_files.split(",")

        if args.main_chain:
            pids_xray = pids_main_chain
        else:
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

            if args.rand_noise:
                # Add noise to the work reflections.
                for i in range(len(f_obs_array.data())):
                    f_ob = f_obs_array.data()[i]
                    f_ob_err = np.random.normal(loc=f_ob, scale=f_ob*.05)
                    f_obs_array.data()[i] = f_ob_err

            # Dropout some work reflections.
            if args.dropout:
                # Define your probabilities
                probabilities = [0.8, 0.2]
                mask = np.random.choice([True, False], size=len(f_obs_array.data()), p=probabilities)
                f_obs_array = f_obs_array.select(flex.bool(mask))

            f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)
            print("N_OBS: ", len(f_obs_array.data()), len(flags_array.data()))

            r_xray = xray_restraint.XtalRestraint(
                hs=hs,
                pids=pids_xray,
                w=w,
                f_obs=f_obs_array,
                free_flags=flags_array,
                w_xray=args.w_xray,
                u_aniso_file=args.u_aniso_file
            )

            rs_xray = IMP.RestraintSet(m, 1.0)
            rs_xray.add_restraint(r_xray)

            r_xrays.append(r_xray)
            rs.append(rs_xray)

            # Setup the weight optimizer state.
            if use_weights and n_states > 1:
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
    step_tracker.set_xray_only(False)
    all_trackers.append(step_tracker)

    time_tracker = trackers.TimeTracker(
        name="time",
        m=m
    )
    time_tracker.set_xray_only(False)
    all_trackers.append(time_tracker)

    ff_tracker = trackers.fTracker(
        name="ff",
        r=rset_charmm
    )
    ff_tracker.set_xray_only(False)
    all_trackers.append(ff_tracker)

    # Add the trackers for each xtal restraint.
    for i in range(len(r_xrays)):
        xray_tracker = trackers.fTracker(
            name="xray_{}".format(i),
            r=r_xrays[i]
        )
        xray_tracker.set_xray_only(True)
        all_trackers.append(xray_tracker)

        r_work_tracker = trackers.RFactorTracker(
            name="r_work_{}".format(i),
            r_xray=r_xrays[i],
            stat="r_work"
        )
        r_work_tracker.set_xray_only(True)
        all_trackers.append(r_work_tracker)

        r_free_tracker = trackers.RFactorTracker(
            name="r_free_{}".format(i),
            r_xray=r_xrays[i],
            stat="r_free"
        )
        r_free_tracker.set_xray_only(True)
        all_trackers.append(r_free_tracker)

    for i in range(len(all_ref_hs)):
        ref_hs = all_ref_hs[i]
        rmsd_all_tracker = trackers.RMSDTracker(
            name="rmsd_avg_{}".format(i),
            hs=hs,
            ref_hs=ref_hs,
            rmsd_func=align_imp.compute_rmsd_between_average,
            ca_only=True
        )
        rmsd_all_tracker.set_xray_only(False)
        all_trackers.append(rmsd_all_tracker)

    for i in range(len(ws)):
        weight_tracker = trackers.WeightTracker(
            name="weight_{}".format(i),
            m=m,
            w=ws[i]
        )
        weight_tracker.set_xray_only(False)
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
    pdb_tracker.set_xray_only(False)
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