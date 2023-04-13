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
sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import molecular_dynamics


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--cif_file")
    parser.add_argument("--res", type=float)
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--dyn_w_xray", type=int)
    parser.add_argument("--restrain_com", type=int)
    parser.add_argument("--start_pdb_file")
    parser.add_argument("--ref_pdb_file")
    parser.add_argument("--T", type=float)
    parser.add_argument("--steps", type=int)
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
    print(args.restrain_com)
    print(args.start_pdb_file)
    print(args.ref_pdb_file)
    print(args.T)
    print(args.steps)
    print(args.sa)
    print(args.log_file)
    print(args.test)

    out_dir = Path(args.out_dir)

    if args.tmp_out_dir == "none":
        tmp_out_dir = None
    else:
        tmp_out_dir = Path(args.tmp_out_dir)

    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=False)

    if tmp_out_dir:
        tmp_out_dir.mkdir(parents=False)

    ref_pdb_file = Path(args.ref_pdb_file)
    pdb_file = Path(args.start_pdb_file)

    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg_symbol = "C 1 2 1"
    cif_file = Path(args.cif_file)

    # if args.test:
    #     steps = 1
    #     pdb_write_freq = None
    # else:
    #     steps = args.steps
    #     pdb_write_freq = 10

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    pids_ca = list(IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())

    if args.res < 0:
        res = None
    else:
        res = args.res

    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmm.charmm_restraints(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)

    r_xtal = xray_restraint.XtalRestraint(
        m=m,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        d_min=res,
        d_max=None,
        scale=True,
        target="ml",
        w_xray=args.w_xray,
        dynamic_w=args.dyn_w_xray
    )

    rs_xtal = IMP.RestraintSet(m, 1.0)
    rs_xtal.add_restraint(r_xtal)

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)

    pids_ca_0 = list(IMP.atom.Selection(h_0, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    com_0 = IMP.atom.CenterOfMass.setup_particle(
        IMP.Particle(m_0),
        pids_ca_0
    )
    print(com_0.get_coordinates())

    rs = [rset_charmm, r_xtal]
    if args.restrain_com:
        r_com = com_restraint.CenterOfMassRestraint(
            m=m,
            pids=pids_ca,
            k=10,
            xyz_0=com_0.get_coordinates()
        )

        rs.append(r_com)

    # rs = [rset_charmm]

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

    rmsd_tracker = trackers.RMSDTracker(
        name="rmsd",
        h=h,
        h_0=h_0,
        align=False
    )
    all_trackers.append(rmsd_tracker)

    rmsd_align_tracker = trackers.RMSDTracker(
        name="rmsd_align",
        h=h,
        h_0=h_0,
        align=True
    )
    all_trackers.append(rmsd_align_tracker)

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

    if args.sa:
        sa_sched = [(300,0,500), (3000,4,50)]
    else:
        sa_sched = None

    log_df = molecular_dynamics.molecular_dynamics(
        output_dir=out_dir,
        tmp_output_dir=tmp_out_dir,
        h=h,
        rs=rs,
        T=args.T,
        t_step=2,
        steps=-1,
        sa_sched=sa_sched,
        all_trackers=all_trackers,
        pdb_write_freq=10,
        log_file=args.log_file
    )