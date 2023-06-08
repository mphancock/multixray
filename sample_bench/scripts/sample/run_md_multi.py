from pathlib import Path
import sys
import argparse
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
import pdb_optimizer_state
import com_optimizer_state
import copy_optimizer_state
import update_weights_optimizer_state
import molecular_dynamics


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--cif_files")
    parser.add_argument("--res", type=float)
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--dyn_w_xray", type=int)
    parser.add_argument("--com")
    parser.add_argument("--start_pdb_file")
    parser.add_argument("--n_state")
    parser.add_argument("--weights", type=int)
    parser.add_argument("--ref_pdb_file")
    parser.add_argument("--T", type=float)
    parser.add_argument("--sa", type=int)
    parser.add_argument("--steps", type=int)
    parser.add_argument("--log_file")
    parser.add_argument("--save_best", action="store_true")
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()
    print("out_dir: ", args.out_dir)
    print("tmp_out_dir: ", args.tmp_out_dir)
    print("cif_files: ", args.cif_files)
    print("res: ", args.res)
    print("w_xray: ", args.w_xray)
    print("dyn_w_xray: ", args.dyn_w_xray)
    print("com: ", args.com)
    print("start_pdb_file: ", args.start_pdb_file)
    print("n_state: ", args.n_state)
    print("weights: ", args.weights)
    print("ref_pdb_file: ", args.ref_pdb_file)
    print("T: ", args.T)
    print("sa: ", args.sa)
    print("steps: ", args.steps)
    print("log_file: ", args.log_file)
    print("save_best: ", args.save_best)
    print("test: ", args.test)

    # Representation.
    pdb_file = Path(args.start_pdb_file)
    ref_pdb_file = Path(args.ref_pdb_file)

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    n_states = int(args.n_state)
    hs = list()
    for i in range(n_states):
        hs.append(IMP.atom.read_pdb(str(pdb_file), m, s))

    print("weights")
    w_p = IMP.Particle(m, "weights")
    w_pid = IMP.isd.Weight.setup_particle(w_p, IMP.algebra.VectorKD([1]*n_states))
    w = IMP.isd.Weight(m, w_pid)

    for i in range(n_states):
        print(w.get_weight(i))

    # Setup waters (if any).
    for h in hs:
        h20s = IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("HET: O  ")).get_selected_particle_indexes()
        for pid in h20s:
            IMP.atom.CHARMMAtom.setup_particle(m, pid, "O")
            print(IMP.atom.CHARMMAtom(m, pid).get_charmm_type())

    ## SAMPLE
    T = args.T
    if args.sa:
        sa_sched = [(T, 0, 500), (3000, 4, 100)]
    else:
        sa_sched = None

    log_file = args.log_file

    pids = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

    save_best = args.save_best

    # We need to manually set the occupancies of the structures here because pdb file occupancies are limited to 2 decimal places.
    for pid in pids:
        IMP.atom.Atom(m, pid).set_occupancy(1/n_states)
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

    r_xrays = list()
    if args.cif_files:
        cif_files = args.cif_files.split(",")
        d_min = args.res
        dyn_w_xray = args.dyn_w_xray
        w_xray = args.w_xray

        # cif file here is a string.
        for cif_file in cif_files:
            r_xray = xray_restraint.XtalRestraint(
                m=m,
                n_state=n_states,
                pids=pids,
                f_obs_file=cif_file,
                d_min=d_min,
                d_max=None,
                scale=True,
                target="ml",
                w_xray=w_xray,
                dynamic_w=dyn_w_xray
            )

            rs_xray = IMP.RestraintSet(m, 1.0)
            rs_xray.add_restraint(r_xray)

            r_xrays.append(r_xray)
            rs.append(rs_xray)

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(Path(ref_pdb_file)), m_0, s_0)

    pids_ca_0 = list(IMP.atom.Selection(h_0, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    com_0 = IMP.atom.CenterOfMass.setup_particle(
        IMP.Particle(m_0),
        pids_ca_0
    )
    print(com_0.get_coordinates())

    all_trackers = list()
    pids_0 = list(IMP.atom.Selection(h_0).get_selected_particle_indexes())

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

    rmsd_tracker = trackers.RMSDTracker(
        name="rmsd",
        hs=hs,
        hs_0=[h_0],
        align=False
    )
    all_trackers.append(rmsd_tracker)

    for i in range(n_states):
        state_pids = IMP.atom.Selection(hs[i]).get_selected_particle_indexes()
        occ_tracker = trackers.OccTracker(
            name="occ_{}".format(i),
            m=m,
            at=IMP.atom.Atom(m, state_pids[0])
        )
        all_trackers.append(occ_tracker)

    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=log_file,
        log_freq=100
    )
    for r in rs:
        r.evaluate(calc_derivs=True)
    log_ostate.update()

    o_states = list()
    o_states.append(log_ostate)

    ## WRITING PDBS
    pdb_dir = Path(args.out_dir, "pdbs")
    pdb_dir.mkdir()

    if args.tmp_out_dir:
        work_pdb_dir = Path(args.tmp_out_dir, "pdbs")
        copy_o_state = copy_optimizer_state.CopyOptimizerState(
            m=m,
            source_dir=work_pdb_dir,
            dest_dir=pdb_dir
        )
        copy_o_state.set_period(100)
        o_states.append(copy_o_state)
    else:
        work_pdb_dir = pdb_dir

    work_pdb_dir.mkdir(exist_ok=True)

    if save_best:
        pdb_ostate = pdb_optimizer_state.WriteBestMultiStatePDBOptimizerState(
            m=m,
            hs=hs,
            pdb_dir=work_pdb_dir,
            tracker=r_free_tracker,
            N=10,
            n_skip=10
        )
    else:
        pdb_ostate = pdb_optimizer_state.WriteMultiStatePDBOptimizerState(
            m=m,
            hs=hs,
            pdb_dir=work_pdb_dir
        )
        pdb_ostate.set_period(10)

    # Write the original frame.
    pdb_ostate.do_update(None)
    o_states.append(pdb_ostate)

    if args.com == "os":
        for h in hs:
            pids_ca = list(IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
            com_os = com_optimizer_state.CenterOfMassOptimizerState(
                m=m,
                pids=pids_ca,
                com_0=com_0
            )
            o_states.append(com_os)

    if args.weights:
        w_os = update_weights_optimizer_state.UpdateWeightsOptimizerState(
            m=m,
            hs=hs,
            w=w,
            r_xtal=r_xtal
        )
        w_os.set_period(10)
        o_states.append(w_os)

    if args.steps:
        steps = args.steps
    else:
        steps = -1

    molecular_dynamics.molecular_dynamics(
        output_dir=work_pdb_dir,
        hs=hs,
        rs=rs,
        T=T,
        t_step=2,
        steps=steps,
        sa_sched=sa_sched,
        o_states=o_states
    )