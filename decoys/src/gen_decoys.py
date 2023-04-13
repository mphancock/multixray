import sys
from pathlib import Path
import shutil
import os

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import log_statistics
import xray_restraint
import charmm
import derivatives
import params


def gen_decoys(
        output_dir,
        h,
        h_0,
        T,
        t_step,
        n_step,
        pdb_write_freq
):
    params.write_params(
        param_dict=locals(),
        param_file=Path(output_dir, "params.txt")
    )

    m = h.get_model()

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    ps = [m.get_particle(pid) for pid in pids]

    rset_charmm = IMP.RestraintSet(m, 1.0)
    rset_charmm.add_restraints(
        charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
    )
    # print(len(m.get_particle_indexes()))
    # for pid in m.get_particle_indexes():
    #     print(m.get_particle_name(pid))

    sf = IMP.core.RestraintsScoringFunction([rset_charmm])

    md = IMP.atom.MolecularDynamics(m)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, T)
    md.add_optimizer_state(s_v)
    md.set_particles(ps)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    # Add all optimizer states.
    trackers = list()
    ff_tracker = log_statistics.fTracker(
        name="ff",
        r=rset_charmm
    )
    trackers.append(ff_tracker)

    m_0 = h_0.get_model()
    pids_0 = IMP.atom.Selection(hierarchy=h_0).get_selected_particle_indexes()
    rmsd_tracker = log_statistics.RMSDTracker(
        name="rmsd",
        m_0=m_0,
        pids_0=pids_0,
        align=True
    )

    log_ostate = log_statistics.LogStatistics(
        m=m,
        pids=pids,
        trackers=[ff_tracker, rmsd_tracker]
    )
    md.add_optimizer_state(log_ostate)

    pdb_dir = Path(output_dir, "pdbs")
    pdb_dir.mkdir(exist_ok=True)
    out_file = Path(pdb_dir, "%1%.pdb")

    if pdb_write_freq:
        pdb_ostate = IMP.atom.WritePDBOptimizerState(
            m,
            pids,
            str(out_file)
        )
        pdb_ostate.set_period(pdb_write_freq)
        md.add_optimizer_state(pdb_ostate)

    md.setup(ps)
    md.set_temperature(T)

    for pid in pids:
        IMP.atom.LinearVelocity.setup_particle(m, pid)

    md.assign_velocities(T)
    md.set_maximum_time_step(t_step)
    md.simulate(n_step * t_step)

    log_df = log_ostate.get_log_as_csv()
    log_df.to_csv(Path(output_dir, "log.csv"))


