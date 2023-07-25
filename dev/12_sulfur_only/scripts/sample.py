import sys
from pathlib import Path
import shutil
import argparse

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import params
import trackers
import log_statistics
import pdb_optimizer_state
import molecular_dynamics
import xray_restraint
import copy_pdbs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--cif_file")
    parser.add_argument("--res", type=float)
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--dyn_w_xray", type=int)
    parser.add_argument("--start_pdb_file")
    parser.add_argument("--ref_pdb_file")
    parser.add_argument("--dof")
    parser.add_argument("--T", type=float)
    parser.add_argument("--sa")
    parser.add_argument("--steps", type=int)
    parser.add_argument("--log_file")
    args = parser.parse_args()
    print("out_dir: ", args.out_dir)
    print("tmp_out_dir: ", args.tmp_out_dir)
    print("cif_file: ", args.cif_file)
    print("res: ", args.res)
    print("w_xray: ", args.w_xray)
    print("dyn_w_xray: ", args.dyn_w_xray)
    print("start_pdb_file: ", args.start_pdb_file)
    print("ref_pdb_file: ", args.ref_pdb_file)
    print("dof: ", args.dof)
    print("T: ", args.T)
    print("sa: ", args.sa)
    print("steps: ", args.steps)
    print("log_file: ", args.log_file)

    pdb_file = Path(args.start_pdb_file)
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    ps = [m.get_particle(pid) for pid in pids]

    out_dir = Path(args.out_dir)
    shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.mkdir(exist_ok=True)
    pdb_dir = Path(out_dir, "pdbs")
    pdb_dir.mkdir()

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

    # XRAY
    if args.cif_file:
        d_min = args.res
        dyn_w_xray = args.dyn_w_xray
        w_xray = args.w_xray

        # cif file here is a string.
        r_xray = xray_restraint.XtalRestraint(
            m=m,
            n_state=1,
            pids=pids,
            f_obs_file=Path(args.cif_file),
            d_min=d_min,
            d_max=None,
            scale=True,
            target="ml",
            w_xray=w_xray,
            dynamic_w=dyn_w_xray
        )

        rs_xray = IMP.RestraintSet(m, 1.0)
        rs_xray.add_restraint(r_xray)
        rs.append(rs_xray)

        xray_tracker = trackers.fTracker(
            name="xray_0",
            r=r_xray
        )
        all_trackers.append(xray_tracker)

        r_work_tracker = trackers.RFactorTracker(
            name="r_work_0",
            r_xray=r_xray,
            stat="r_work"
        )
        all_trackers.append(r_work_tracker)

        r_free_tracker = trackers.RFactorTracker(
            name="r_free_0",
            r_xray=r_xray,
            stat="r_free"
        )
        all_trackers.append(r_free_tracker)

    ff_tracker = trackers.fTracker(
        name="ff",
        r=rset_charmm
    )
    all_trackers.append(ff_tracker)

    m_0 = IMP.Model()
    h_0 = IMP.atom.read_pdb(str(args.ref_pdb_file), m_0, IMP.atom.AllPDBSelector())
    rmsd_tracker = trackers.RMSDTracker(
        name="rmsd",
        hs=[h],
        hs_0=[h_0],
        align=False
    )
    all_trackers.append(rmsd_tracker)

    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=str(Path(out_dir, "log.csv")),
        log_freq=10
    )

    # pdb_ostate = pdb_optimizer_state.WriteMultiStatePDBOptimizerState(
    #         m=m,
    #         hs=[h],
    #         pdb_dir=pdb_dir
    # )
    # pdb_ostate.set_period(10)
    # pdb_ostate.do_update(None)

    pids = list()
    pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())
    # ps = [m.get_particle(pid) for pid in pids]
    ps = list()


    # for pid in pids:
    #     IMP.core.XYZ(m, pid).set_coordinates_are_optimized(False)

    # for pid in pids_235:
    #     IMP.core.XYZ(m, pid).set_coordinates_are_optimized(True)

    # cg = IMP.core.ConjugateGradients(m)
    # cg.set_scoring_function(sf)
    # cg.add_optimizer_state(log_ostate)
    # cg.add_optimizer_state(pdb_ostate)
    # cg.optimize(100)

    # for pid in pids:
    #     IMP.atom.LinearVelocity.setup_particle(m, pid)

    # md = IMP.atom.MolecularDynamics(m)
    # md.set_particles(ps_235)
    # md.set_scoring_function(sf)
    # md.set_has_required_score_states(True)

    # T=4000
    # md.setup(ps)
    # md.set_temperature(T)
    # s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, T)
    # md.add_optimizer_state(s_v)
    # md.add_optimizer_state(log_ostate)
    # md.add_optimizer_state(pdb_ostate)

    # md.assign_velocities(T)
    # md.set_maximum_time_step(2)

    # md.simulate(1e99)
    if args.sa:
        sa_sched = args.sa.split(";")
        sa_sched = [[float(T), float(res), int(steps), dof_str] for T, res, steps, dof_str in [x.split(",") for x in sa_sched]]
        for i in range(len(sa_sched)):
            float, res, steps, dof_str = sa_sched[i]

            if dof_str == "A":
                pids_work = pids
            elif dof_str == "235_side":
                pids_work = (IMP.atom.Selection(hierarchy=h, residue_index=235) - IMP.atom.Selection(hierarchy=h, atom_types=[IMP.atom.AT_CA, IMP.atom.AT_C, IMP.atom.AT_O, IMP.atom.AT_N])).get_selected_particle_indexes()
            elif dof_str == "235":
                pids_work = (IMP.atom.Selection(hierarchy=h, residue_index=235)).get_selected_particle_indexes()
            elif dof_str == "S_side":
                pids_work = (IMP.atom.Selection(hierarchy=h, residue_types=[IMP.atom.MET, IMP.atom.CYS]) - IMP.atom.Selection(hierarchy=h, atom_types=[IMP.atom.AT_CA, IMP.atom.AT_C, IMP.atom.AT_O, IMP.atom.AT_N])).get_selected_particle_indexes()
            elif dof_str == "S":
                pids_work = (IMP.atom.Selection(hierarchy=h, residue_types=[IMP.atom.MET, IMP.atom.CYS])).get_selected_particle_indexes()
            else:
                raise RuntimeError()

            sa_sched[i][3] = pids_work
    else:
        sa_sched = None

    ps_work = [m.get_particle(pid) for pid in pids_work]
    ps = [m.get_particle(pid) for pid in pids]

    if args.tmp_out_dir:
        work_pdb_dir = Path(args.tmp_out_dir, "pdbs")
        copy_o_state = copy_pdbs.CopyOptimizerState(
            m=m,
            source_dir=work_pdb_dir,
            dest_dir=pdb_dir
        )
        copy_o_state.set_period(100)
    else:
        work_pdb_dir = pdb_dir

    work_pdb_dir.mkdir(exist_ok=True)

    pdb_ostate = pdb_optimizer_state.WriteMultiStatePDBOptimizerState(
        m=m,
        hs=[h],
        pdb_dir=work_pdb_dir
    )
    pdb_ostate.set_period(10)
    pdb_ostate.do_update(None)

    molecular_dynamics.molecular_dynamics(
        output_dir=out_dir,
        hs=[h],
        rs=rs,
        T=args.T,
        t_step=2,
        steps=-1,
        sa_sched=sa_sched,
        o_states=[log_ostate, pdb_ostate, copy_o_state],
        md_ps=ps
    )
