import time

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

    def get_name(
            self
    ):
        return self.name

    def get_model(
            self
    ):
        return self.m

    def get_n(
            self
    ):
        return self.n

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
        return tuple(self.xyz.get_coordinates())


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
        return tuple(self.xyz.get_derivatives())


class RMSDTracker(Tracker):
    def __init__(
            self,
            name,
            hs,
            hs_0,
            align
    ):
        Tracker.__init__(
            self,
            name=name,
            m=hs[0].get_model(),
            n=1
        )
        self.hs = hs
        self.hs_0 = hs_0
        self.align = align

    def evaluate(
            self
    ):
        rmsd = align_imp.compute_rmsd_between_average(
            hs_1=self.hs,
            hs_2=self.hs_0
        )

        return rmsd


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
        # Try to see if the restraint has a get_f call that stores the score.
        try:
            return self.r.get_f()
        except AttributeError:
            return self.r.evaluate(False)


# Stat must be one of the following: r_free, r_work, r_all.
class RFactorTracker(Tracker):
    def __init__(
            self,
            name,
            r_xray,
            stat
    ):
        Tracker.__init__(
            self,
            name=name,
            m=r_xray.get_model(),
            n=1
        )
        self.r_xray = r_xray
        self.stat = stat

    def evaluate(
            self
    ):
        if self.stat == "r_free":
            return self.r_xray.get_r_free()
        elif self.stat == "r_work":
            return self.r_xray.get_r_work()
        elif self.stat == "r_all":
            return self.r_xray.get_r_all()
        else:
            raise RuntimeError("Incorrect stat type. Must be one of the following: r_free, r_work, r_all.")


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
        return t


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

    def evaluate(
            self
    ):
        step = self.step
        self.step = self.step + 1
        return step