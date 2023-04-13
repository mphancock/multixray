import IMP
import IMP.core
import IMP.atom


class CenterOfMassOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            pids,
            com_0
    ):
        IMP.OptimizerState.__init__(self, m, "CenterOfMassOptimizerState%1%")
        self.m = m
        self.pids = pids
        self.com_0 = com_0

    def do_update(self, call):
        com = IMP.atom.CenterOfMass.setup_particle(
            IMP.Particle(self.m),
            self.pids
        )

        delta = self.com_0.get_coordinates() - com.get_coordinates()
        for pid in self.pids:
            xyz = IMP.core.XYZ(self.m, pid)
            coords_adjust = xyz.get_coordinates() + delta
            xyz.set_coordinates(coords_adjust)


