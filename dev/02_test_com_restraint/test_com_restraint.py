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

    out_dir = Path(Path.home(), "xray/dev/02_test_com_restraint/output_0")

    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=False)

    ref_pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    pdb_file = ref_pdb_file
    # pdb_file = Path(Path.home(), "xray/dev/02_test_com_restraint/3ca7_clean_drift.pdb")

    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg_symbol = "C 1 2 1"

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    pids_ca = list(IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())

    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmm.charmm_restraints(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)
    pids_ca_0 = list(IMP.atom.Selection(h_0, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    com_0 = IMP.atom.CenterOfMass.setup_particle(
        IMP.Particle(m_0),
        pids_ca_0
    )
    print(com_0.get_coordinates())

    r_com = com_restraint.CenterOfMassRestraint(
        m=m,
        pids=pids_ca,
        k=10,
        xyz_0=com_0.get_coordinates()
    )

    cif_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")
    r_xtal = xray_restraint.XtalRestraint(
        m=m,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        d_min=4,
        d_max=None,
        scale=True,
        target="ml",
        w_xray=30000,
        dynamic_w=0
    )

    # rs = [rset_charmm, r_xtal]
    rs = [rset_charmm, r_com, r_xtal]

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

    log_df = molecular_dynamics.molecular_dynamics(
        output_dir=out_dir,
        tmp_output_dir=None,
        h=h,
        rs=rs,
        T=3000,
        t_step=2,
        steps=5000,
        sa_sched=None,
        all_trackers=all_trackers,
        pdb_write_freq=10,
        log_file=None
    )