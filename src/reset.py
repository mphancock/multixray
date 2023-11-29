import numpy as np
from pathlib import Path

import IMP
import IMP.core
import IMP.atom

import trackers


"""
Doesn't work perfectly. The system has a tendency to quickly jump back to high energy states following the rest.
"""
class ResetTracker(trackers.Tracker):
    def __init__(
            self,
            name,
            hs,
            w,
            r_charmm,
            sa_sched,
            h_0
    ):
        trackers.Tracker.__init__(
            self,
            name=name,
            m=hs[0].get_model(),
            n=1
        )

        self.hs = hs
        self.w = w
        self.r_charmm = r_charmm
        self.md = None
        self.sv = None
        self.sa_sched = sa_sched
        self.h_0 = h_0
        self.step = 0

    def get_md(self):
        return self.md

    def set_md(
            self,
            md
    ):
        self.md = md

    def set_sv(
            self,
            sv
    ):
        self.sv = sv

    def reset_positions(self):
        pids_0 = IMP.atom.Selection(self.h_0).get_selected_particle_indexes()
        for h in self.hs:
            pids = IMP.atom.Selection(h).get_selected_particle_indexes()

            if len(pids_0) != len(pids):
                raise RuntimeError("The number of particles in the initial state and the current state are not the same.")

            for i in range(len(pids)):
                xyz = IMP.core.XYZ(self.m, pids[i])
                xyz_0 = IMP.core.XYZ(self.h_0.get_model(), pids_0[i])

                xyz.set_coordinates(xyz_0.get_coordinates())

                # Set the derivatives and velocities to 0.
                da = IMP.DerivativeAccumulator()
                d = IMP.core.XYZR(self.get_model(), pids[0])
                d.add_to_derivatives(-d.get_derivatives(), da)

                IMP.atom.LinearVelocity(self.m, pids[i]).set_velocity([0,0,0])


    def determine_current_T(self):
        current_step = self.step % np.sum([self.sa_sched[i]["step"] for i in range(len(self.sa_sched))])

        for i in range(len(self.sa_sched)):
            steps = self.sa_sched[i]["step"]

            if current_step < steps:
                return self.sa_sched[i]["T"]
            else:
                current_step = current_step - steps

        raise RuntimeError("Something has gone wrong")


    def reset_velocities(self):
        T = self.determine_current_T()
        T=300

        # Need to take care of the scaling optimizer and the md engine.
        print("RESETTING TO {}K".format(T))
        self.md.set_temperature(T)
        self.sv.set_temperature(T)
        # self.md.assign_velocities(T)


    def do_evaluate(self):
        self.reset_positions()
        self.reset_velocities()

        return 1

    def evaluate(
            self
    ):
        # ff_score = self.r_charmm.evaluate(False)
        ff_score = self.r_charmm.get_last_score()

        # We need this first condition because the tracker will be evaluated before the MD has been set.
        if self.step > 0 and self.step % self.get_period() == 0 and self.get_on() and ff_score > 1e5 and self.get_md():
            reset = self.do_evaluate()
            IMP.atom.write_multimodel_pdb(self.hs, str(Path(Path.home(), "xray/tmp/tmp.pdb")))
            # ff_score = self.r_charmm.evaluate(False)
        else:
            reset = 0

        if self.step > 0:
            pids = IMP.atom.Selection(self.hs[0]).get_selected_particle_indexes()
            print("T: {}".format(self.md.get_temperature()))
            print("KE: {}".format(self.md.get_kinetic_energy()))
            print("ff: {}".format(ff_score))
            print("gradients: {}".format(IMP.core.XYZR(self.m, pids[0]).get_derivatives()))
            print("velocity: {}".format(IMP.atom.LinearVelocity(self.m, pids[0]).get_velocity()))

        self.step = self.step+1

        return reset
