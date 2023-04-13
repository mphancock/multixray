from pathlib import Path
import sys
import argparse
import shutil

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import xray_restraint
import log_statistics
import derivatives
sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import molecular_dynamics
import merge_pdbs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("--job_name")
    # parser.add_argument("--job_id", type=int)
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--")
    args = parser.parse_args()
    # print(args.job_name)
    # print(args.job_id)
    # # print(args.cif_file)
    # # print(args.w_xray)
    print(args.test)

    job_name = "single_sample_volume"
    job_id = 0

    job_dir = Path(Path.home(), "xray/sample_bench/out/3ca7/{}".format(job_name))
    output_dir = Path(job_dir, "output_{}".format(job_id))

    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=False)

    ref_pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    pdb_file = Path(Path.home(), "xray/decoys/data/decoy_sets/3ca7/3ca7_N_1000_x1/315.pdb")
    # pdb_file = ref_pdb_file

    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg_symbol = "C 1 2 1"
    w_xray = 10**4
    cif_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")

    steps = 2000
    T = 3000
    pdb_write_freq = 4

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    r_xtal = xray_restraint.XtalRestraint(
        m=m,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        d_min=4,
        d_max=None,
        scale=True,
        target="ml"
    )
    r_xtal.set_weight(w_xray)

    rs_xtal = IMP.RestraintSet(m, 1.0)
    rs_xtal.add_restraint(r_xtal)
    # rs_xtal.set_weight(w_xray)

    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmm.charmm_restraints(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)

    trackers = list()
    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)
    pids_0 = list(IMP.atom.Selection(h_0).get_selected_particle_indexes())

    rmsd_tracker = log_statistics.RMSDTracker(
        name="RMSD",
        h_0=h_0,
        align=True
    )
    trackers.append(rmsd_tracker)

    ff_tracker = log_statistics.fTracker(
        name="ff",
        r=rset_charmm
    )
    trackers.append(ff_tracker)

    xray_tracker = log_statistics.fTracker(
        name="xray",
        r=r_xtal
    )
    trackers.append(xray_tracker)

    df_mag_tracker = log_statistics.dfMagTracker(
        name="df_mag",
        r1=rset_charmm,
        r2=r_xtal
    )
    trackers.append(df_mag_tracker)

    dxray_tracker = log_statistics.dfTracker(
        name="AT1 dxraydx",
        r=r_xtal,
        pid=pids[0],
        at_id=0
    )
    trackers.append(dxray_tracker)

    molecular_dynamics.molecular_dynamics(
        output_dir=output_dir,
        h=h,
        rs=[r_xtal, rset_charmm],
        T=T,
        t_step=2,
        steps=steps,
        trackers=trackers,
        pdb_write_freq=pdb_write_freq,
        test=args.test
    )

    merge_pdbs(
        output_dir=output_dir,
        prefix_pdb_file=ref_pdb_file
    )