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


def sim_anneal(
        output_dir,
        h,
        uc_dim,
        sg_symbol,
        cif_file,
        w_xray,
        t_step,
        n_cycles,
        n_steps,
        temp_schedule,
        d_min_schedule,
        pdb_freq,
        h_0
):
    if len(temp_schedule) != len(d_min_schedule):
        raise RuntimeError("Simulated annealing temperature {} and resolution {} schedules must be of equal size!".format(len(temp_schedule), len(d_min_schedule)))

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
    s_v = IMP.atom.VelocityScalingOptimizerState(
        m=m,
        pis=ps,
        temperature=temp_schedule[0]
    )
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

    temp_tracker = log_statistics.tempTracker(
        name="T",
        md=md
    )

    log_ostate = log_statistics.LogStatistics(
        m=m,
        pids=pids,
        trackers=[temp_tracker, ff_tracker, xray_tracker, df_mag_tracker, rmsd_tracker]
    )
    md.add_optimizer_state(log_ostate)

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

    md.set_time_step(t_step)

    for i in range(n_cycles):
        for j in range(len(temp_schedule)):
            T = temp_schedule[j]
            d_min = d_min_schedule[j]

            r_xtal.set_d_min(
                d_min=d_min
            )
            s_v.set_temperature(T)
            md.set_temperature(T)
            # md.assign_velocities(T)
            md.set_maximum_time_step(t_step)

            md.simulate(n_steps * t_step)

    # SIMULATION
    log_df = log_ostate.get_log_as_csv()
    log_df.to_csv(Path(output_dir, "log.csv"))