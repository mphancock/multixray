import sys
from pathlib import Path

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import params
import reset
import pdb_writer


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
"""
def molecular_dynamics(
        pdb_dir,
        hs,
        r_sets,
        t_step,
        n_step,
        sa_sched,
        o_states
):
    params.write_params(
        param_dict=locals(),
        param_file=Path(pdb_dir.parents[0], "params_md.txt")
    )

    m = hs[0].get_model()

    pids = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

    ps_0 = [m.get_particle(pid) for pid in sa_sched[0]["dof"]]
    T_0 = sa_sched[0]["T"]

    md = IMP.atom.MolecularDynamics(m)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps_0, T_0)
    md.add_optimizer_state(s_v)
    md.set_particles(ps_0)
    sf = IMP.core.RestraintsScoringFunction(r_sets)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    for o_state in o_states:
        md.add_optimizer_state(o_state)

    md.setup(ps_0)
    md.set_temperature(T_0)
    for pid in pids:
        IMP.atom.LinearVelocity.setup_particle(m, pid)

    md.assign_velocities(T_0)
    md.set_maximum_time_step(t_step)

    # Set MD in the reset tracker.
    trackers = o_states[0].get_trackers()
    for tracker in trackers:
        if type(tracker) == reset.ResetTracker:
            tracker.set_md(md)

    # rs is a set of restraint sets
    if sa_sched:
        # Get the log ostate, always exists.
        log_ostate = o_states[0]

        write_pdb_tracker = None
        copy_pdb_tracker = None
        for tracker in trackers:
            print(tracker)
            if type(tracker) == pdb_writer.PDBWriterTracker:
                write_pdb_tracker = tracker
            elif type(tracker) == pdb_writer.PDBCopyTracker:
                copy_pdb_tracker = tracker

        # Get the r_xrays, if there is one.
        r_names = [r_sets[i].get_restraint(0).get_name() for i in range(len(r_sets))]
        r_set_xrays = list()
        r_set_not_xrays = list()
        for i in range(len(r_names)):
            r = r_sets[i].get_restraint(0)
            if "XrayRestraint" in r_names[i]:
                # Add the restraint set not the restraint.
                r_set_xrays.append(r_sets[i])
            else:
                r_set_not_xrays.append(r_sets[i])

        # Get the scoring function that exludes the xray restraints.
        sf = IMP.core.RestraintsScoringFunction(r_sets)
        sf_not_xray = IMP.core.RestraintsScoringFunction(r_set_not_xrays)

        # Get the weight ostate, if exists.
        weight_ostate = None
        for o_state in o_states:
            print(o_state.get_name())
            if "UpdateWeightsOptimizerState" in o_state.get_name():
                weight_ostate = o_state

        step_tracker = log_ostate.get_tracker("step")
        while step_tracker.get_step() < n_step:
            for sa_step in sa_sched:
                print(sa_step)

                n_step_sa = sa_step["step"]
                T = sa_step["T"]
                pids_dof = sa_step["dof"]
                pdb = sa_step["pdb"]
                w = sa_step["w"]
                res = sa_step["res"]

                md.set_particles([m.get_particle(pid) for pid in pids_dof])

                # If d_min is negative, then turn off the xray restraint.
                if res < 0:
                    md.set_scoring_function(sf_not_xray)

                    # Turn off xray only.
                    for tracker in log_ostate.get_trackers():
                        if tracker.get_xray_only():
                            print("Turning off: {}".format(tracker.get_name()))
                            tracker.set_off()
                else:
                    md.set_scoring_function(sf)

                    for r_set_xray in r_set_xrays:
                        r_xray = r_set_xray.get_restraint(0)
                        r_xray.set_d_min(
                            d_min=res
                        )

                    # Turn on xray only.
                    for tracker in log_ostate.get_trackers():
                        if tracker.get_xray_only():
                            print("Turning on: {}".format(tracker.get_name()))
                            tracker.set_on()

                if w:
                    if weight_ostate:
                        weight_ostate.set_on(True)
                else:
                    if weight_ostate:
                        weight_ostate.set_on(False)

                if pdb:
                    write_pdb_tracker.set_on()
                    if copy_pdb_tracker:
                        copy_pdb_tracker.set_on()
                else:
                    write_pdb_tracker.set_off()
                    if copy_pdb_tracker:
                        copy_pdb_tracker.set_off()

                s_v.set_temperature(T)
                md.set_temperature(T)
                md.set_maximum_time_step(t_step)

                md.simulate(n_step_sa*t_step)
    else:
        md.simulate(n_step*t_step)

    # IMP.atom.write_multimodel_pdb(hs, str(Path(output_dir, "%1.pdb")))
    return 0
