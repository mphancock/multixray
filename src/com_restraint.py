

import IMP
import IMP.core
import IMP.atom


"""
The CenterOfMassRestraint is a restraint on the center of mass of the model. The restraint is implemented as a simple harmonic with spring coefficient k. The restraint is useful for keeping the model from drifting. 
"""
class CenterOfMassRestraint(IMP.Restraint):
    def __init__(
            self,
            m,
            pids,
            k,
            xyz_0
    ):
        IMP.Restraint.__init__(
            self,
            m,
            "CenterOfMassRestraint%1%"
        )
        # The model.
        self.m = m
        # All atomic model pids.
        self.pids = pids
        # The spring coefficient.
        self.k = k
        # The coordinates of the reference model center of mass.
        self.xyz_0 = xyz_0

    def do_add_score_and_derivatives(self, sa):
        com = IMP.atom.CenterOfMass.setup_particle(
            IMP.Particle(self.m),
            self.pids
        )

        delta = (com.get_coordinates()-self.xyz_0)
        mag_delta = delta.get_magnitude()
        score = .5*mag_delta**2

        if sa.get_derivative_accumulator():
            deriv = self.k*delta
            # com.add_to_derivative(
            #     deriv,
            #     sa.get_derivative_accumulator()
            # )

            for pid in self.pids:
                d = IMP.core.XYZR(self.m, pid)
                d.add_to_derivatives(
                    deriv,
                    sa.get_derivative_accumulator()
                )

        sa.add_score(score)

    def do_get_inputs(self):
        return [self.get_model().get_particle(pid) for pid in self.pids]