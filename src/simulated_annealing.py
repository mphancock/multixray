import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import params
import reset
import pdb_writer
import log_statistics


class SimulatedAnnealingSchedule:
    def __init__(
        self,
        sa_string
    ):
        self.read_sa_sched_string(sa_string)
        print(self.sa_sched_df)

    def read_sa_sched_string(
        self,
        sa_str
    ):
        self.sa_sched_df = pd.DataFrame

        sa_sched_strs = sa_str.split(";")

        sa_sched = defaultdict(list)

        ## s is number of steps, t is temperature, d is degrees of freedom, p is pdb, w is weight, r is resolution
        keys = ["s", "t", "d", "p", "w", "r"]
        for sa_step_str in sa_sched_strs:
            for key_val_pair in sa_step_str.split(","):
                key = key_val_pair[0]
                val = key_val_pair[1:]

                if key in ["s", "t", "p", "w"]:
                    val = int(val)
                elif key == "r":
                    val = float(val)

                sa_sched[key].append(val)

        self.sa_sched_df = pd.DataFrame(sa_sched)

    def get_n_sa_steps(
        self
    ):
        return len(self.sa_sched_df)

    def get_steps(
        self,
        sa_step
    ):
        return self.sa_sched_df.loc[sa_step, "s"]

    def get_temperature(
        self,
        sa_step
    ):
        return self.sa_sched_df.loc[sa_step, "t"]

    def get_degrees_of_freedom(
        self,
        sa_step
    ):
        return self.sa_sched_df.loc[sa_step, "d"]

    def get_pdb(
        self,
        sa_step
    ):
        return self.sa_sched_df.loc[sa_step, "p"]

    def get_weight(
        self,
        sa_step
    ):
        return self.sa_sched_df.loc[sa_step, "w"]

    def get_resolution(
        self,
        sa_step
    ):
        return self.sa_sched_df.loc[sa_step, "r"]

    def get_weight_on(
        self,
        sa_step
    ):
        return self.get_weight(sa_step) == 1

    def get_pdb_on(
        self,
        sa_step
    ):
        return self.get_pdb(sa_step) == 1

    def get_xray_on(
        self,
        sa_step
    ):
        return self.get_resolution(sa_step) >= 0

    def get_pids_for_degrees_of_freedom(
        self,
        sa_step,
        msmc_m
    ):
        dof = self.get_degrees_of_freedom(sa_step)

        if dof == "A":
            pids_dof = msmc_m.get_all_pids()
        elif dof == "S":
            pids_dof = msmc_m.get_all_side_pids()
        elif dof[0].isnumeric():
            start, end = dof.split("-")
            pids_dof = msmc_m.get_pids_in_res_range(int(start), int(end))
        else:
            raise ValueError("Invalid dof: {}".format(dof))

        return pids_dof


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
class SimulatedAnnealing:
    def __init__(
        self,
        msmc_m,
        rset_xray,
        rset_charmm,
        t_step,
        n_step,
        sa_sched,
        log_o_state,
        weight_o_state,
        com_o_state
    ):
        self.msmc_m = msmc_m
        self.rset_xray = rset_xray
        self.rset_charmm = rset_charmm
        self.t_step = t_step
        self.n_step = n_step
        self.sa_sched = sa_sched
        self.log_o_state = log_o_state
        self.weight_o_state = weight_o_state
        self.com_o_state = com_o_state

        self.hs = msmc_m.get_hs()
        self.m = msmc_m.get_m()

        self.pids = msmc_m.get_all_pids()
        self.ps = [self.m.get_particle(pid) for pid in self.pids]

        # Setup the md
        ## Reorder in a bid...
        self.md = IMP.atom.MolecularDynamics(self.m)
        T_0 = sa_sched.get_temperature(0)
        self.s_v = IMP.atom.VelocityScalingOptimizerState(self.m, self.ps, T_0)
        self.md.add_optimizer_state(self.s_v)
        self.md.set_particles(self.ps)

        sf = IMP.core.RestraintsScoringFunction([rset_charmm])
        self.md.set_scoring_function(sf)
        self.md.set_has_required_score_states(True)

        for o_state in [self.log_o_state, self.weight_o_state, self.com_o_state]:
            if o_state:
                self.md.add_optimizer_state(o_state)

        self.md.setup(self.ps)
        self.md.set_temperature(T_0)
        for pid in self.pids:
            IMP.atom.LinearVelocity.setup_particle(self.m, pid)

        self.md.assign_velocities(T_0)
        self.md.set_maximum_time_step(self.t_step)

    def run(self):
        cur_step = 0
        cur_sa_step = 0
        while cur_step < self.n_step:
            ## turn on/off the xray restraints
            if self.sa_sched.get_xray_on(cur_sa_step):
                ## set the resolution of all xray restraints
                res = self.sa_sched.get_resolution(cur_sa_step)
                for i in range(self.rset_xray.get_number_of_restraints()):
                    r_xray = self.rset_xray.get_restraint(i)
                    print(type(r_xray))
                    r_xray.set_d_min(res)
                sf = IMP.core.RestraintsScoringFunction([self.rset_xray, self.rset_charmm])

                ## turn on the xray related trackers (eg, r free)
                for tracker in self.log_o_state.get_trackers():
                    if tracker.get_xray_only():
                        tracker.set_on()
            else:
                sf = IMP.core.RestraintsScoringFunction([rset_charmm])

                ## turn off the xray related trackers (eg, r free)
                for tracker in self.log_o_state.get_trackers():
                    if tracker.get_xray_only():
                        tracker.set_off()

            self.md.set_scoring_function(sf)

            ## set the pids for degrees of freedom
            self.md.set_particles([self.m.get_particle(pid) for pid in self.sa_sched.get_pids_for_degrees_of_freedom(cur_sa_step, self.msmc_m)])

            ## turn on/off the weight optimizer states
            self.weight_o_state.set_on(self.sa_sched.get_weight_on(cur_sa_step))

            ## turn on/off the pdb writing
            write_pdb_tracker = self.log_o_state.get_trackers_by_type(pdb_writer.PDBWriterTracker)[0]
            if self.sa_sched.get_pdb_on(cur_sa_step):
                write_pdb_tracker.set_on()
            else:
                write_pdb_tracker.set_off()

            T = self.sa_sched.get_temperature(cur_sa_step)
            n_frames = self.sa_sched.get_steps(cur_sa_step)
            t_step = 2

            self.s_v.set_temperature(T)
            self.md.set_temperature(T)
            self.md.set_maximum_time_step(t_step)
            # md.set_velocity_cap(.005)

            self.md.simulate(n_frames*t_step)

            cur_step = cur_step + n_frames
            cur_sa_step = (cur_sa_step + 1) % sa_sched.get_n_sa_steps()

        return 0
