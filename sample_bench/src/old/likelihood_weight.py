from pathlib import Path
import sys
import os
import shutil

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
sys.path.append(str(Path(Path.home(), "xray/src")))
import log_statistics
import xray_restraint
import charmm
import params


def likelihood_weight(
        output_dir,
        h,
        uc_dim,
        sg_symbol,
        cif_file,
        w_xray,
        t_step,
        steps,
        T,
        equil,
        pdb_freq,
        h_0
):
    params.write_params(
        param_dict=locals(),
        param_file=Path(output_dir, "params.txt")
    )

    m = h.get_model()
    m_0 = h_0.get_model()

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    ps = [m.get_particle(pid) for pid in pids]
    r_xtal = xray_restraint.XtalRestraint(
        m=m,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        d_max=None,
        d_min=4,
        scale=True,
        target="ml"
    )

    rs_xtal = IMP.RestraintSet(m, 1.0)
    rs_xtal.add_restraint(r_xtal)

    rset_charmm = IMP.RestraintSet(m, 1.0)
    rset_charmm.add_restraints(
        charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
    )

    rs_xtal.set_weight(w_xray)
    sf = IMP.core.RestraintsScoringFunction([rs_xtal, rset_charmm])

    md = IMP.atom.MolecularDynamics(m)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, T)
    md.add_optimizer_state(s_v)
    md.set_particles(ps)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    for pid in pids:
        IMP.atom.LinearVelocity.setup_particle(m, pid)

    md.setup(ps)

    # Add all optimizer states.
    pids_track = IMP.atom.Selection(
        hierarchy=h,
        atom_type=IMP.atom.AtomType("CA")
    ).get_selected_particle_indexes()

    # m_0 = IMP.Model()
    # s_0 = IMP.atom.AllPDBSelector()
    # h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)
    pids_0 = list(IMP.atom.Selection(h_0).get_selected_particle_indexes())

    ff_tracker = log_statistics.fTracker(
        name="ff",
        r=rset_charmm
    )

    xray_tracker = log_statistics.fTracker(
        name="xray",
        r=r_xtal
    )

    rmsd_tracker = log_statistics.RMSDTracker(
        name="RMSD",
        m_0=m_0,
        pids_0=pids_0
    )

    df_mag_tracker = log_statistics.dfMagTracker(
        name="df_mag",
        r1=rset_charmm,
        r2=r_xtal
    )

    m2_sd_s = IMP.atom.Selection(
        hierarchy=h,
        residue_index=2,
        atom_type=IMP.atom.AtomType("SD")
    )
    m2_sd_pid = m2_sd_s.get_selected_particle_indexes()[0]

    dff_tracker = log_statistics.dfTracker(
        name="M2 SD dff dx",
        r=rset_charmm,
        pid=m2_sd_pid,
        at_id=0
    )

    dxray_tracker = log_statistics.dfTracker(
        name="M2 SD dxray dx",
        r=r_xtal,
        pid=pids[0],
        at_id=0
    )

    log_ostate = log_statistics.LogStatistics(
        m=m,
        pids=pids,
        trackers=[ff_tracker, xray_tracker, df_mag_tracker, rmsd_tracker, dff_tracker, dxray_tracker]
    )
    md.add_optimizer_state(log_ostate)

    # WARMUP
    # if equil:
    # md.set_temperature(100)
    # md.assign_velocities(100)
    # md.set_maximum_time_step(2)
    # md.simulate(50)

    pdb_dir = Path(output_dir, "pdbs")
    pdb_dir.mkdir(exist_ok=True)
    out_file = Path(pdb_dir, "%1%.pdb")
    pdb_ostate = IMP.atom.WritePDBOptimizerState(
        m,
        pids,
        str(out_file)
    )
    pdb_ostate.set_period(pdb_freq)
    md.add_optimizer_state(pdb_ostate)

    # SIMULATION
    md.set_temperature(T)
    md.assign_velocities(T)

    lv = IMP.atom.LinearVelocity(m, pid)
    print(lv.get_velocity())

    vs_os = IMP.atom.VelocityScalingOptimizerState(
        m=m,
        pis=pids,
        temperature=T
    )
    vs_os.set_period(10)
    md.add_optimizer_state(vs_os)

    md.set_maximum_time_step(t_step)
    md.simulate(steps * t_step)

    log_df = log_ostate.get_log_as_csv()
    log_df.to_csv(Path(output_dir, "log.csv"))


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")
    output_dir = Path(sys.argv[1])

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    Path(output_dir).mkdir(parents=False)
    Path(output_dir, "pdbs").mkdir(parents=False)

    ref_pdb_file = Path(Path.home(), "xray/data/pdbs/4n7f/4n7f_heavy.pdb")

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(sys.argv[4]), m, s)

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)

    cif_file = Path(sys.argv[3])
    sg_symbol = "P 41 2 2"
    uc_dim = (68.411, 68.411, 37.248, 90.0, 90.0, 90.0)

    steps = 2500
    T = 5000
    t_step = 2
    w_xray = 1000
    pdb_freq = 5

    if Path(sys.argv[4]).name != "4n7f_heavy.pdb":
        equil = False
    else:
        equil = True

    likelihood_weight(
        output_dir=output_dir,
        h=h,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        cif_file=cif_file,
        w_xray=w_xray,
        t_step=t_step,
        steps=steps,
        T=T,
        equil=equil,
        pdb_freq=pdb_freq,
        h_0=h_0
    )

