import sys
from pathlib import Path

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import log_statistics
import copy_optimizer_state
import params
import pdb_writers


def molecular_dynamics(
        output_dir,
        hs,
        rs,
        T,
        t_step,
        steps,
        sa_sched,
        o_states
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
    ps = [m.get_particle(pid) for pid in pids]

    sf = IMP.core.RestraintsScoringFunction(rs)

    md = IMP.atom.MolecularDynamics(m)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, T)
    md.add_optimizer_state(s_v)
    md.set_particles(ps)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    pdb_dir = Path(output_dir, "pdbs")

    for o_state in o_states:
        md.add_optimizer_state(o_state)

    md.setup(ps)
    md.set_temperature(T)
    for pid in pids:
        IMP.atom.LinearVelocity.setup_particle(m, pid)

    md.assign_velocities(T)
    md.set_maximum_time_step(t_step)

    # sf.evaluate(derivatives=True)
    # log_ostate.update()

    if sa_sched:
        while True:
            for T, d_min, steps in sa_sched:
                rs[1].set_d_min(
                    d_min=d_min
                )
                s_v.set_temperature(T)
                md.set_temperature(T)
                # md.assign_velocities(T)
                md.set_maximum_time_step(t_step)

                md.simulate(steps*t_step)
    else:
        if steps > 0:
            md.simulate(steps*t_step)
        else:
            md.simulate(1e99)

    # log_df = log_ostate.get_log()

    return 0
