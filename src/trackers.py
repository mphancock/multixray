import time
from typing import Any
import numpy as np
import shutil
import os

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
            n,
            labels=None
    ):
        self.name = name
        self.m = m
        self.n = n
        self.on = True
        self.xray_only = False
        self.period = 1

        if labels:
            self.labels = labels
        else:
            if self.n > 1:
                self.labels = ["{}_{}".format(name, i) for i in range(self.n)]
            else:
                self.labels = [name]

        if len(self.labels) != self.n:
            raise ValueError("Number of labels ({}) does not match n ({})".format(len(self.labels), self.n))

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
            hs_0,
            hs_1,
            pids_0,
            pids_1,
            occs_0,
            occs_1
    ):
        Tracker.__init__(
            self,
            name=name,
            m=hs_0[0].get_model(),
            n=1
        )
        self.hs_0 = hs_0
        self.hs_1 = hs_1
        self.pids_0 = pids_0
        self.pids_1 = pids_1
        self.occs_0 = occs_0
        self.occs_1 = occs_1
        self.rmsd_func = rmsd_func


    def evaluate(
            self
    ):
        rmsd = self.rmsd_func(
            h_0s=self.hs_0,
            h_1s=self.hs_1,
            pids_0=self.pids_0,
            pids_1=self.pids_1,
            occs_0=self.occs_0,
            occs_1=self.occs_1
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
            r_xray,
            labels
    ):
        Tracker.__init__(
            self,
            name=name,
            m=r_xray.get_model(),
            n=2,
            labels=labels
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
            pid,
            scale

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
        self.scale = scale

    def evaluate(
            self
    ):
        df_dict = derivatives.get_df_dict(
            m=self.m,
            pids=self.pids,
            r=self.r
        )
        return df_dict[self.pid] * self.scale


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


class WeightMatTracker(Tracker):
    def __init__(
            self,
            name,
            msmc_m,
            labels
    ):
        Tracker.__init__(
            self,
            name=name,
            m=msmc_m.get_m(),
            n=msmc_m.get_w_mat().shape[0]*msmc_m.get_w_mat().shape[1],
            labels=labels
        )
        self.msmc_m = msmc_m

        n_weights = msmc_m.get_w_mat().shape[0] * msmc_m.get_w_mat().shape[1]

        self.set_labels(labels)

    def evaluate(self):
        return list(self.msmc_m.get_w_mat().flatten())


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


class CopyTracker(Tracker):
    def __init__(
            self,
            name,
            m,
            source_dir,
            dest_dir
    ):
        Tracker.__init__(
            self,
            name=name,
            m=m,
            n=1
        )
        self.source_dir = source_dir
        self.dest_dir = dest_dir
        self.step = 0

    def do_evaluate(self):
        # for file in self.source_dir.glob("*.pdb"):
        #     dest_file = Path(self.dest_dir, file.name)
        #     shutil.move(file, dest_file)

        # shutil.copytree(self.source_dir, self.dest_dir)

        # Iterate over the source directory contents
        print("Performing copy after {} steps".format(self.step))
        for item in os.listdir(self.source_dir):
            s = os.path.join(self.source_dir, item)
            d = os.path.join(self.dest_dir, item)
            if os.path.isdir(s):
                # Recursively copy subdirectory
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                # Copy file
                shutil.copy2(s, d)

        return 1

    def evaluate(
            self
    ):
        if self.step % self.get_period() == 0 and self.get_on():
            copied = self.do_evaluate()
        else:
            copied = 0

        self.step = self.step+1

        return [copied]