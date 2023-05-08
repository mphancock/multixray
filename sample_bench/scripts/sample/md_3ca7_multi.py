from pathlib import Path
import sys
import argparse
import shutil

import IMP
import IMP.atom
import IMP.core

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import xray_restraint
import trackers
import com_restraint
import log_statistics
import align_imp
import pdb_writers
import com_optimizer_state
import copy_optimizer_state
import molecular_dynamics


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--cif_file")
    parser.add_argument("--res", type=float)
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--dyn_w_xray", type=int)
    parser.add_argument("--com")
    parser.add_argument("--start_pdb_file")
    parser.add_argument("--ref_pdb_file")
    parser.add_argument("--T", type=float)
    parser.add_argument("--sa", type=int)
    parser.add_argument("--log_file")
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()
    print(args.out_dir)
    print(args.tmp_out_dir)
    print(args.cif_file)
    print(args.res)
    print(args.w_xray)
    print(args.dyn_w_xray)
    print(args.com)
    print(args.start_pdb_file)
    print(args.ref_pdb_file)
    print(args.T)
    print(args.sa)
    print(args.log_file)
    print(args.test)

    output_dir = args.out_dir
    if args.tmp_out_dir == "none":
        tmp_out_dir = None
    else:
        tmp_out_dir = Path(args.tmp_out_dir)

    pdb_file = Path(args.start_pdb_file)
    ref_pdb_file = Path(args.ref_pdb_file)

    # Representation.
    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, s)
    n_states = len(hs)

    # Scoring.
    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg_symbol = "C 1 2 1"
    cif_file = Path(args.cif_file)
    d_min = args.res
    dyn_w_xray = args.dyn_w_xray
    w_xray = args.w_xray

    # Sampling.
    T = args.T
    sa = args.sa

    log_file = args.log_file

    pids = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

    # We need to manually set the occupancies of the structures her because pdb file occupancies are limited to 2 decimal places.
    for pid in pids:
        IMP.atom.Atom(m, pid).set_occupancy(1/n_states)
    ps = [m.get_particle(pid) for pid in pids]

    rset_charmm = IMP.RestraintSet(m, 1.0)
    for h in hs:
        charmm_rs = charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
        rset_charmm.add_restraints(charmm_rs)
    rset_charmm.set_weight(1/n_states)

    r_xtal = xray_restraint.XtalRestraint(
        m=m,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        d_min=d_min,
        d_max=None,
        scale=True,
        target="ml",
        w_xray=w_xray,
        dynamic_w=dyn_w_xray
    )

    rs_xtal = IMP.RestraintSet(m, 1.0)
    rs_xtal.add_restraint(r_xtal)

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")), m_0, s_0)

    pids_ca_0 = list(IMP.atom.Selection(h_0, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    com_0 = IMP.atom.CenterOfMass.setup_particle(
        IMP.Particle(m_0),
        pids_ca_0
    )
    print(com_0.get_coordinates())

    tmp_pdb_dir = Path(tmp_out_dir, "pdbs")
    pdb_dir = Path(output_dir, "pdbs")
    tmp_pdb_dir.mkdir()
    pdb_dir.mkdir()

    o_states = list()
    pdb_ostate = pdb_writers.WriteMultiStatePDBOptimizerState(
        m=m,
        hs=hs,
        pdb_dir=tmp_pdb_dir
    )
    pdb_ostate.update()
    pdb_ostate.set_period(10)
    o_states.append(pdb_ostate)

    copy_o_state = copy_optimizer_state.CopyOptimizerState(
        m=m,
        source_dir=tmp_pdb_dir,
        dest_dir=pdb_dir
    )
    copy_o_state.set_period(100)
    o_states.append(copy_o_state)

    if args.com == "os":
        for h in hs:
            pids_ca = list(IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
            com_os = com_optimizer_state.CenterOfMassOptimizerState(
                m=m,
                pids=pids_ca,
                com_0=com_0
            )
            o_states.append(com_os)

    rs = [rset_charmm, r_xtal]
    # for h in hs:
    #     pids_ca = list(IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    #
    #     r_com = com_restraint.CenterOfMassRestraint(
    #         m=m,
    #         pids=pids_ca,
    #         k=10,
    #         xyz_0=com_0.get_coordinates()
    #     )
    #
    #     rs.append(r_com)

    # sf = IMP.core.RestraintsScoringFunction(rs)

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

    xray_tracker = trackers.fTracker(
        name="xray",
        r=r_xtal
    )
    all_trackers.append(xray_tracker)

    r_work_tracker = trackers.RFactorTracker(
        name="r_work",
        r_xray=r_xtal,
        stat="r_work"
    )
    all_trackers.append(r_work_tracker)

    r_free_tracker = trackers.RFactorTracker(
        name="r_free",
        r_xray=r_xtal,
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

    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=log_file,
        log_freq=100
    )
    o_states.append(log_ostate)

    for r in rs:
        r.evaluate(calc_derivs=True)
    log_ostate.update()

    molecular_dynamics.molecular_dynamics(
        output_dir=tmp_out_dir,
        hs=hs,
        rs=rs,
        T=T,
        t_step=2,
        steps=-1,
        sa_sched=sa,
        o_states=o_states
    )