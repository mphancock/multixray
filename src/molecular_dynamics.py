import sys
from pathlib import Path

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import params


"""
Run a generic multi-state molecular dynamics simulation.

**********
Parameters:
    output_dir: the output directory for the simulation containing the params file, log file, and pdb files.

    hs: a list of hierarchies to simulate.

    rs: a list of restraint sets to use for scoring.

    T: the temperature to use for the simulation.

    t_step: the time step to use for the simulation.

    steps: the number of steps to run the simulation for. A T=-1 means run until the process is killed.

    sa_sched: a list of tuples of the form (T, d_min, steps, pids_work) where T is the temperature to use for the simulation, d_min is the resolution cutoff to use for the xray restraint, steps is the number of steps to run the simulation for, and pids_work is the list of particle ids to use for the simulation.

    o_states: a list of optimizer states to use for the simulation.

    md_ps: a list of particles to use for the simulation. If None, then all particles in the hierarchies are used. This is redundant with the pids_work argument in sa_sched, but is included for convenience.
"""
def molecular_dynamics(
        output_dir,
        hs,
        rs,
        T,
        t_step,
        steps,
        sa_sched,
        o_states,
        md_ps=None
):
    params.write_params(
        param_dict=locals(),
        param_file=Path(output_dir, "params.txt")
    )

    m = hs[0].get_model()
    n_states = len(hs)
    print("N_STATES: {}".format(n_states))

    pids = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

    if md_ps:
        ps = md_ps
    else:
        ps = [m.get_particle(pid) for pid in pids]

    sf = IMP.core.RestraintsScoringFunction(rs)

    md = IMP.atom.MolecularDynamics(m)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, T)
    md.add_optimizer_state(s_v)
    md.set_particles(ps)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    for o_state in o_states:
        md.add_optimizer_state(o_state)

    md.setup(ps)
    md.set_temperature(T)
    for pid in pids:
        IMP.atom.LinearVelocity.setup_particle(m, pid)

    md.assign_velocities(T)
    md.set_maximum_time_step(t_step)

    # rs is a set of restraint sets
    if sa_sched:
        r_xray = rs[1].get_restraint(0)

        log_ostate = o_states[0]
        xray_only_tracker_names = ["xray_0", "r_work_0", "r_free_0", "pdb", "copy"]
        xray_only_trackers = list()
        for name in xray_only_tracker_names:
            xray_only_trackers.append(log_ostate.get_tracker(name))

        while True:
            for T, d_min, steps, pids_work in sa_sched:
                print(T, d_min, steps)

                ps_work = [m.get_particle(pid) for pid in pids_work]
                md.set_particles(ps_work)

                if d_min < 0:
                    # Need to turn off dynamic xray scaling as well.
                    sf = IMP.core.RestraintsScoringFunction([rs[0]])
                    md.set_scoring_function(sf)

                    # Turn off xray and pdb writing.
                    for tracker in xray_only_trackers:
                        tracker.set_writing(False)
                else:
                    r_xray.set_d_min(
                        d_min=d_min
                    )
                    sf = IMP.core.RestraintsScoringFunction(rs)
                    md.set_scoring_function(sf)

                    for tracker in xray_only_trackers:
                        tracker.set_writing(True)

                s_v.set_temperature(T)
                md.set_temperature(T)
                md.set_maximum_time_step(t_step)

                md.simulate(steps*t_step)
    else:
        if steps > 0:
            md.simulate(steps*t_step)
        else:
            md.simulate(1e99)

    # log_df = log_ostate.get_log()

    return 0
