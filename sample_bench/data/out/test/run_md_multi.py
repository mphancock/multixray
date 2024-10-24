from pathlib import Path
import sys
import argparse
import shutil
import random

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
import copy_pdbs
import update_weights_optimizer_state
import molecular_dynamics
import params
import weight_restraint
import miller_ops


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--cif_files")
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--dyn_w_xray", action="store_true")
    parser.add_argument("--start_pdb_file")
    parser.add_argument("--n_state")
    parser.add_argument("--init_weights")
    parser.add_argument("--weights", action="store_true")
    parser.add_argument("--ref_pdb_file")
    parser.add_argument("--T", type=float)
    parser.add_argument("--sa")
    parser.add_argument("--steps", type=int)
    parser.add_argument("--save_best", action="store_true")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--bfactor", type=int)
    parser.add_argument("--dropout", action="store_true")
    parser.add_argument("--d_min", type=float)
    parser.add_argument("--rand_noise", action="store_true")
    parser.add_argument("--main_chain", action="store_true")
    args = parser.parse_args()

    params.write_params(vars(args), Path(args.out_dir, "params.txt"))

    # Representation.
    pdb_file = Path(args.start_pdb_file)
    ref_pdb_file = Path(args.ref_pdb_file)

    m = IMP.Model()

    s = IMP.atom.AllPDBSelector()

    n_states = int(args.n_state)
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, s)

    if len(hs) == 1 and len(hs) < n_states:
        m = IMP.Model()
        hs.clear()

        for i in range(n_states):
            hs.append(IMP.atom.read_pdb(str(pdb_file), m, s))

    pids = list()

    m_0 = IMP.Model()
    h0_s = IMP.atom.read_multimodel_pdb(str(ref_pdb_file), m_0, s)

    # Get all particle ids.
    pids = list()
    pids_main_chain = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())
        pids_main_chain.extend(IMP.atom.Selection(h, atom_types=[IMP.atom.AT_CA, IMP.atom.AT_C, IMP.atom.AT_O, IMP.atom.AT_N]).get_selected_particle_indexes())

    # Setup bfactors.
    if args.bfactor:
        for h in hs:
            for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
                IMP.atom.Atom(m, pid).set_temperature_factor(args.bfactor)

    # Setup the weights.
    w_p = IMP.Particle(m, "weights")
    w_pid = IMP.isd.Weight.setup_particle(w_p, IMP.algebra.VectorKD([1]*n_states))
    w = IMP.isd.Weight(m, w_pid)

    if args.init_weights == "ref":
        init_weights = [IMP.atom.Atom(m, IMP.atom.Selection(hierarchy=h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()[0]).get_occupancy() for h in hs]
    elif args.init_weights == "rand":
        init_weights = [random.random() for i in range(n_states)]
    elif args.init_weights == "uni":
        init_weights = [1/n_states]*n_states
    else:
        init_weights = args.init_weights.split(",")
        init_weights = [float(w) for w in init_weights]
        if len(init_weights) != n_states:
            raise RuntimeError("Number of initial weights does not match number of states.")

    # Normalize and check for trivial weights only if init weights was rand. I do this because some of the synthetic natives have trivial weights but I don't want to overrule them.
    if args.init_weights == "rand":
        weights = [w/sum(init_weights) for w in init_weights]
        trivial_w = False
        for weight in weights:
            if weight < .05:
                trivial_w = True
        while trivial_w:
            weights = [random.random() for i in range(n_states)]
            weights = [w/sum(weights) for w in weights]
            trivial_w = False
            for weight in weights:
                if weight < .05:
                    trivial_w = True
    else:
        weights = init_weights

    w.set_weights(weights)

    # We need to manually set the occupancies of the structures here because pdb file occupancies are limited to 2 decimal places.
    for i in range(n_states):
        pids = IMP.atom.Selection(hs[i]).get_selected_particle_indexes()
        for pid in pids:
            IMP.atom.Atom(m, pid).set_occupancy(w.get_weight(i))

    w.set_weights_are_optimized(True)

    # Setup waters (if any).
    for h in hs:
        h20s = IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("HET: O  ")).get_selected_particle_indexes()
        for pid in h20s:
            IMP.atom.CHARMMAtom.setup_particle(m, pid, "O")

    ## SAMPLE
    T = args.T

    save_best = args.save_best

    if args.sa:
        sa_sched = args.sa.split(";")
        sa_sched = [[float(T), float(res), int(steps), dof_str] for T, res, steps, dof_str in [x.split(",") for x in sa_sched]]
        for i in range(len(sa_sched)):
            float, res, steps, dof_str = sa_sched[i]

            if dof_str == "A":
                pids_work = pids
            elif dof_str == "S_side":
                pids_work = (IMP.atom.Selection(hierarchy=h, residue_types=[IMP.atom.MET, IMP.atom.CYS]) - IMP.atom.Selection(hierarchy=h, atom_types=[IMP.atom.AT_CA, IMP.atom.AT_C, IMP.atom.AT_O, IMP.atom.AT_N])).get_selected_particle_indexes()
            elif dof_str == "S":
                pids_work = (IMP.atom.Selection(hierarchy=h, residue_types=[IMP.atom.MET, IMP.atom.CYS])).get_selected_particle_indexes()
            else:
                raise RuntimeError()

            sa_sched[i][3] = pids_work
    else:
        sa_sched = None

    ps = [m.get_particle(pid) for pid in pids]

    ## SCORING
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
    if args.cif_files:
        cif_files = args.cif_files.split(",")
        if args.main_chain:
            pids_xray = pids_main_chain
        else:
            pids_xray = pids

        # cif file here is a string.
        for cif_file in cif_files:
            r_xray = xray_restraint.XtalRestraint(
                m=m,
                n_state=n_states,
                pids=pids_xray,
                f_obs_file=cif_file,
                scale=True,
                target="ml",
                w_xray=args.w_xray,
                dynamic_w=args.dyn_w_xray,
                dropout=args.dropout,
                d_min=args.d_min,
                rand_noise=args.rand_noise
            )

            rs_xray = IMP.RestraintSet(m, 1.0)
            rs_xray.add_restraint(r_xray)

            r_xrays.append(r_xray)
            rs.append(rs_xray)

            if args.weights:
                f_obs = miller_ops.get_miller_array(
                    f_obs_file=cif_file,
                    label="_refln.F_meas_au"
                )
                f_obs_array = miller_ops.clean_miller_array(f_obs)
                # Set flags.
                status_array = miller_ops.get_miller_array(
                    f_obs_file=cif_file,
                    label="_refln.status"
                )
                flags_array = status_array.customized_copy(data=status_array.data()=="f")
                f_obs, flags_array = f_obs_array.common_sets(other=flags_array)

                wr = weight_restraint.WeightRestraint(
                    m=m,
                    hs=hs,
                    pids=pids_xray,
                    w=w,
                    f_obs=f_obs,
                    flags=flags_array,
                    scale=0.5
                )
                wrs.append(wr)

    # Setup the center of mass optimizer state.
    pids_ca_0 = list(IMP.atom.Selection(h0_s, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    com_0 = IMP.atom.CenterOfMass.setup_particle(
        IMP.Particle(m_0),
        pids_ca_0
    )

    all_trackers = list()
    pids_0 = list(IMP.atom.Selection(h0_s).get_selected_particle_indexes())

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
        all_trackers.append(xray_tracker)

        r_work_tracker = trackers.RFactorTracker(
            name="r_work_{}".format(i),
            r_xray=r_xrays[i],
            stat="r_work"
        )
        all_trackers.append(r_work_tracker)

        r_free_tracker = trackers.RFactorTracker(
            name="r_free_{}".format(i),
            r_xray=r_xrays[i],
            stat="r_free"
        )
        all_trackers.append(r_free_tracker)

    rmsd_all_tracker = trackers.RMSDTracker(
        name="rmsd_avg",
        hs=hs,
        hs_0=h0_s,
        rmsd_func=align_imp.compute_rmsd_between_average,
        ca_only=True
    )
    all_trackers.append(rmsd_all_tracker)

    for i in range(n_states):
        state_pids = IMP.atom.Selection(hs[i]).get_selected_particle_indexes()
        occ_tracker = trackers.OccTracker(
            name="occ_{}".format(i),
            m=m,
            at=IMP.atom.Atom(m, state_pids[0])
        )
        all_trackers.append(occ_tracker)

    ## WRITING PDBS
    pdb_dir = Path(args.out_dir, "pdbs")
    pdb_dir.mkdir()

    tmp_pdb_dir = Path(args.tmp_out_dir, "pdbs")
    tmp_pdb_dir.mkdir()

    pdb_tracker = pdb_writer.PDBWriterTracker(
        name="pdb",
        hs=hs,
        pdb_dir=tmp_pdb_dir,
        freq=10,
        log_pdb_dir=pdb_dir
    )
    all_trackers.append(pdb_tracker)

    # This freq needs to be equal to the logging frequency to ensure there are never more log entries than corresponding pdb files.
    copy_tracker = copy_pdbs.PDBCopyTracker(
        name="copy",
        m=m,
        source_dir=tmp_pdb_dir,
        dest_dir=pdb_dir,
        freq=100
    )
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

    o_states = list()
    o_states.append(log_ostate)

    for h in hs:
        pids_ca = list(IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
        com_os = com_optimizer_state.CenterOfMassOptimizerState(
            m=m,
            pids=pids_ca,
            com_0=com_0
        )
        o_states.append(com_os)

    # Setup the weight optimizer state.
    if args.weights and n_states > 1:
        # w_os = update_weights_optimizer_state.OptimizeWeightsOptimizerState(
        #     m=m,
        #     wrs=wrs,
        #     w=w,
        #     n_state=n_states,
        #     step_tracker=step_tracker
        # )

        w_os = update_weights_optimizer_state.UpdateWeightsOptimizerState(
            m=m,
            hs=hs,
            w=w,
            r_xrays=r_xrays,
            n_proposals=10,
            radius=.05
        )

        w_os.set_period(100)
        o_states.append(w_os)

    if args.steps:
        steps = args.steps
    else:
        steps = 1e99

    molecular_dynamics.molecular_dynamics(
        pdb_dir=tmp_pdb_dir,
        hs=hs,
        rs=rs,
        T=T,
        t_step=2,
        steps=steps,
        sa_sched=sa_sched,
        o_states=o_states
    )

    pdb_tracker.evaluate()
    copy_tracker.evaluate()