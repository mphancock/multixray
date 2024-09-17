import IMP
import IMP.core
import IMP.atom


class CenterOfMassOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            pids,
            ref_com
    ):
        IMP.OptimizerState.__init__(self, m, "CenterOfMassOptimizerState%1%")
        self.m = m
        self.pids = pids
        self.ref_com = ref_com

        self.is_on = True

    def get_is_on(self):
        return self.is_on

    def turn_on(self):
        self.is_on = on

    def turn_off(self):
        self.is_on = False

    def do_update(self, call):
        if self.is_on:
            com = IMP.atom.CenterOfMass.setup_particle(
                IMP.Particle(self.m),
                self.pids
            )

            delta = self.ref_com.get_coordinates() - com.get_coordinates()
            for pid in self.pids:
                xyz = IMP.core.XYZ(self.m, pid)
                coords_adjust = xyz.get_coordinates() + delta
                xyz.set_coordinates(coords_adjust)


