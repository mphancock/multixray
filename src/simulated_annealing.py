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
from trackers import TempTracker
from weight_thermostat import WeightThermostat


class SimulatedAnnealingSchedule:
    def __init__(
        self,
        sa_string
    ):
        self.read_sa_sched_string(sa_string)

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
            pids_dof = msmc_m.get_pids()
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
    msmc_m (MultiStateMultiConditionModel): The multi-state multi-condition model.

    rset_xray (IMP.RestraintSet): The x-ray restraints.

    rset_charmm (IMP.RestraintSet): The CHARMM restraints.

    t_step (float): The time step.

    n_step (int): The number of steps.

    sa_sched (SimulatedAnnealingSchedule): The simulated annealing schedule.

    log_o_state (IMP.OptimizerState): The log optimizer state.

    weight_o_state (IMP.OptimizerState): The weight optimizer state.

    com_o_state (IMP.OptimizerState): The center of mass optimizer state.
"""
class SimulatedAnnealing:
    def __init__(
        self,
        msmc_m,
        rset_xray,
        r_charmm,
        t_step,
        n_step,
        sa_sched,
        log_o_state,
        weight_o_state,
        com_o_state,
        vel_thermo,
        weight_thermo
    ):
        self.msmc_m = msmc_m
        self.rset_xray = rset_xray
        self.r_charmm = r_charmm
        self.t_step = t_step
        self.n_step = n_step
        self.sa_sched = sa_sched
        self.log_o_state = log_o_state
        self.weight_o_state = weight_o_state
        self.com_o_state = com_o_state

        self.hs = msmc_m.get_hs()
        self.m = msmc_m.get_m()

        self.pids = msmc_m.get_pids()
        self.ps = [self.m.get_particle(pid) for pid in self.pids]

        # Setup the md
        ## Reorder in a bid...
        self.md = IMP.atom.MolecularDynamics(self.m)
        T_0 = sa_sched.get_temperature(0)

        # self.s_v = IMP.atom.VelocityScalingOptimizerState(self.m, self.ps, T_0)
        # self.s_v.set_period(10)
        self.md.set_particles(self.ps)

        self.ff_sf = IMP.core.RestraintsScoringFunction([r_charmm])
        self.xray_sf = IMP.core.RestraintsScoringFunction([r_charmm, rset_xray])

        self.md.set_scoring_function(self.ff_sf)
        self.md.set_has_required_score_states(True)

        # setup velocity thermostat but only on non solvent atoms
        ## this doesn't work bc the particle stack on the thermostat and md must be the same
        # self.prot_pids = msmc_m.get_all_protein_pids()
        # self.prot_ps = [self.m.get_particle(pid) for pid in self.prot_pids]

        if vel_thermo:
            self.vel_thermo_o_state = IMP.atom.BerendsenThermostatOptimizerState(
                pis=self.ps,
                temperature=T_0,
                tau=10
            )
        else:
            self.vel_thermo_o_state = None

        ## create weight thermostat
        if rset_xray.get_number_of_restraints() > 0:
            self.weight_thermostat_o_state = WeightThermostat(
                msmc_m=msmc_m,
                rset_xray=rset_xray,
                md=self.md,
                T_target=T_0,
                warmup_steps=50
            )
            self.weight_thermostat_o_state.set_period(10)

            if not weight_thermo:
                self.weight_thermostat_o_state.set_on(False)
        else:
            self.weight_thermostat_o_state = None

        for o_state in [self.com_o_state, self.weight_o_state, self.vel_thermo_o_state, self.weight_thermostat_o_state, self.log_o_state]:
            if o_state:
                self.md.add_optimizer_state(o_state)

        ## turn off center of mass adjustment
        self.com_o_state.turn_off()

        self.md.setup(self.ps)
        self.md.set_temperature(T_0)
        self.md.assign_velocities(T_0)

        # set solvent velocities to 0
        # for pid in msmc_m.get_all_solvent_pids():
        #     IMP.atom.LinearVelocity(self.m, pid).set_velocity(IMP.algebra.Vector3D(0,0,0))

        self.md.set_maximum_time_step(self.t_step)

        ## if there is a temp tracker turn it on
        temp_tracker = self.log_o_state.get_trackers_by_type(TempTracker)[0]
        temp_tracker.set_on()
        temp_tracker.set_md(self.md)


    def run(self):
        cur_step = 0
        cur_sa_step = 0
        while cur_step < self.n_step:
            print(self.sa_sched.get_weight_on(cur_sa_step))
            ## turn on/off the xray restraints
            if self.sa_sched.get_xray_on(cur_sa_step):
                ## set the resolution of all xray restraints
                res = self.sa_sched.get_resolution(cur_sa_step)
                for i in range(self.rset_xray.get_number_of_restraints()):
                    r_xray = self.rset_xray.get_restraint(i)
                    r_xray.set_d_min(res)

                ## need to evaluate in this order
                # sf = IMP.core.RestraintsScoringFunction([self.rset_charmm, self.shadow_restraint, self.rset_xray])

                ## turn on the xray related trackers (eg, r free)
                for tracker in self.log_o_state.get_trackers():
                    if tracker.get_xray_only():
                        tracker.set_on()

                self.md.set_scoring_function(self.xray_sf)
            else:
                # sf = IMP.core.RestraintsScoringFunction([rset_charmm])

                ## turn off the xray related trackers (eg, r free)
                for tracker in self.log_o_state.get_trackers():
                    if tracker.get_xray_only():
                        tracker.set_off()

                self.md.set_scoring_function(self.ff_sf)

            # self.md.set_scoring_function(sf)

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

            # self.s_v.set_temperature(T)
            self.md.set_temperature(T)
            self.md.set_maximum_time_step(t_step)
            self.md.set_velocity_cap(.25)

            self.md.simulate(n_frames*t_step)

            cur_step = cur_step + n_frames
            cur_sa_step = (cur_sa_step + 1) % sa_sched.get_n_sa_steps()

        return 0
