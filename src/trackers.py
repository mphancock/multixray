import time
from typing import Any
import numpy as np

import IMP
import IMP.atom
import IMP.core

import derivatives
import align_imp


class Tracker:
    def __init__(
            self,
            name,
            m,
            n
    ):
        self.name = name
        self.m = m
        self.n = n
        self.on = True
        self.xray_only = False
        self.period = 1

        if n > 1:
            self.labels = ["{}_{}".format(name, i) for i in range(self.n)]
        else:
            self.labels = [name]

    def get_name(self):
        return self.name

    def get_model(self):
        return self.m

    def get_n(self):
        return self.n

    def get_on(self):
        return self.on

    def get_xray_only(self):
        return self.xray_only

    def get_period(self):
        return self.period

    def get_labels(self):
        return self.labels

    def set_name(
            self,
            new_name
    ):
        self.name = new_name

    def set_model(
            self,
            model
    ):
        self.m = model

    def set_n(
            self,
            n
    ):
        self.n = n

    def set_on(self):
        self.on = True

    def set_off(self):
        self.on = False

    def set_xray_only(
            self,
            xray_only
    ):
        self.xray_only = xray_only

    def set_period(
            self,
            period
    ):
        self.period = period

    def set_labels(
            self,
            labels
    ):
        self.labels = labels


class XYZTracker(Tracker):
    def __init__(
            self,
            name,
            m,
            xyz
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=3
        )
        self.xyz = xyz

    def evaluate(
            self
    ):
        return list(self.xyz.get_coordinates())


class dXYZTracker(Tracker):
    def __init__(
            self,
            name,
            m,
            xyz
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=3
        )
        self.xyz = xyz

    def evaluate(
            self
    ):
        return list(self.xyz.get_derivatives())


class RMSDTracker(Tracker):
    def __init__(
            self,
            name,
            rmsd_func,
            hs,
            w,
            ref_hs,
            ref_w,
            ca_only
    ):
        Tracker.__init__(
            self,
            name=name,
            m=hs[0].get_model(),
            n=1
        )
        self.hs = hs
        self.hs_0 = ref_hs
        self.w = w
        self.ref_w = ref_w
        self.rmsd_func = rmsd_func
        self.ca_only = ca_only

    def evaluate(
            self
    ):
        rmsd = self.rmsd_func(
            h_0s=self.hs,
            h_1s=self.hs_0,
            w_0=self.w,
            w_1=self.ref_w,
            ca_only=self.ca_only
        )

        return [float(rmsd)]


class fTracker(Tracker):
    def __init__(
            self,
            name,
            r

    ):
        Tracker.__init__(
            self,
            name=name,
            m=r.get_model(),
            n=1
        )
        # self.r can be a restraint or a restraint set.
        self.r = r

    def evaluate(
            self
    ):
        if not self.get_on():
            return np.nan

        # Try to see if the restraint has a get_f call that stores the score.
        try:
            last_score = self.r.get_f()
        except AttributeError:
            last_score = self.r.get_last_score()

        return [last_score]


class RFactorTracker(Tracker):
    def __init__(
            self,
            name,
            r_xray
    ):
        Tracker.__init__(
            self,
            name=name,
            m=r_xray.get_model(),
            n=2
        )
        self.r_xray = r_xray
        self.set_xray_only(True)

    def evaluate(
            self
    ):
        if not self.get_on():
            r_free, r_work = np.nan, np.nan
        else:
            r_free, r_work = self.r_xray.get_r_free(), self.r_xray.get_r_work()

        return [r_free, r_work]


class dfdXYZTracker(Tracker):
    def __init__(
            self,
            name,
            m,
            pids,
            r,
            pid

    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=3
        )
        # self.r can be a restraint or a restraint set.
        self.r = r
        self.pid = pid
        self.pids = pids

    def evaluate(
            self
    ):
        df_dict = derivatives.get_df_dict(
            m=self.m,
            pids=self.pids,
            r=self.r
        )
        return df_dict[self.pid]


class dfMagTracker(Tracker):
    def __init__(
            self,
            name,
            m,
            pids,
            r1,
            r2
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=1
        )
        # r1 and r2 can be a restraint or a restraint set.
        self.r1 = r1
        self.r2 = r2
        self.pids = pids

    def evaluate(
            self
    ):
        ratio = derivatives.get_df_mag_ratio(
            m=self.m,
            pids=self.pids,
            r1=self.r1,
            r2=self.r2
        )

        return ratio


class TempTracker(Tracker):
    def __init__(
            self,
            name,
            md
    ):
        Tracker.__init__(
            self,
            name=name,
            m=md.get_model(),
            n=1
        )
        self.md = md

    def evaluate(
            self
    ):
        # energy in kcal / mol
        energy = self.md.get_kinetic_energy()
        temp = self.md.get_kinetic_temperature(
            ekinetic=energy
        )

        return temp

class OccTracker(Tracker):
    def __init__(
            self,
            name,
            m,
            at
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=1
        )
        self.at = at

    def evaluate(
            self
    ):
        occ = self.at.get_occupancy()

        return occ

class WeightTracker(Tracker):
    def __init__(
            self,
            name,
            m,
            w
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=w.get_number_of_weights()
        )
        self.w = w

    def evaluate(self):
        return list(self.w.get_weights())

class TimeTracker(Tracker):
    def __init__(
            self,
            name,
            m
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=1
        )
        self.t_0 = None

    def evaluate(
            self
    ):
        if not self.t_0:
            self.t_0 = time.time()

        t = time.time() - self.t_0
        return [t]


class StepTracker(Tracker):
    def __init__(
            self,
            name,
            m
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=1
        )
        self.step = 0

    def get_step(
            self
    ):
        return self.step

    def evaluate(
            self
    ):
        step = self.step
        self.step = self.step + 1
        return [step]


class DFMagRatioTracker(Tracker):
    def __init__(
            self,
            name,
            m,
            xr
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=1
        )

        self.xr = xr

    def evaluate(
            self
    ):
        return self.xr.get_df_mag_ratio()