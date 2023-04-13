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
sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import molecular_dynamics
import merge_pdbs


def get_w_xray(
        job_id
):
    if 0 <= job_id <= 9:
        w_xray = (job_id)*100+100
    elif 10 <= job_id <= 29:
        w_xray = (job_id-9)*200+1000
    elif 30 <= job_id <= 69:
        w_xray = (job_id-29)*500+5000
    elif 70 <= job_id <= 79:
        w_xray = (job_id-70)*.1+.1
    elif 80 <= job_id <= 99:
        w_xray = (job_id-79)*.2+1
    elif 100 <= job_id <= 139:
        w_xray = (job_id-99)*.5+5

    return w_xray


def get_dynamic_w(
        job_id
):
    if job_id <= 69:
        return 0
    else:
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_out_dir")
    parser.add_argument("--job_id", type=int)
    parser.add_argument("--run_id", type=int)
    args = parser.parse_args()

    print(args.out_dir)
    print(args.tmp_out_dir)
    print(args.job_id)
    print(args.run_id)

    out_dir = Path(args.out_dir)
    tmp_out_dir = Path(args.tmp_out_dir)
    log_file = Path(out_dir, "log.csv")

    ref_pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    pdb_file = ref_pdb_file

    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg_symbol = "C 1 2 1"
    cif_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")
    res = -1
    steps = 999999999

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()

    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmm.charmm_restraints(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)

    w_xray = get_w_xray(job_id=args.job_id)
    dynamic_w = get_dynamic_w(job_id=args.job_id)

    print("w_xray: {}".format(w_xray))
    print("dynamic_w: {}".format(dynamic_w))

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
        w_xray=w_xray,
        dynamic_w=dynamic_w

    )
    rs_xtal = IMP.RestraintSet(m, 1.0)

    # rs = [rset_charmm, r_xtal]
    rs = [rset_charmm]

    all_trackers = list()
    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)
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
        align=True
    )
    all_trackers.append(rmsd_tracker)

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

    log_df = molecular_dynamics.molecular_dynamics(
        output_dir=out_dir,
        tmp_output_dir=tmp_out_dir,
        h=h,
        rs=rs,
        T=300,
        t_step=2,
        steps=steps,
        all_trackers=all_trackers,
        pdb_write_freq=10,
        log_file=log_file
    )