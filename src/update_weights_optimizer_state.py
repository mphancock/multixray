import IMP
import IMP.core
import IMP.atom


class UpdateWeightsOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            hs,
            w,
            r_xtal
    ):
        IMP.OptimizerState.__init__(self, m, "UpdateWeightsOptimizerState%1%")
        self.m = m
        self.hs = hs
        self.w = w
        self.r_xtal = r_xtal

    def do_update(self, call):
        # First, update w from w_grads.
        w_grads = self.r_xtal.get_w_grads()

        w_grads_scale = [w_grads[i]*1000 for i in range(len(w_grads))]
        ws = self.w.get_weights()
        for i in range(len(ws)):
            ws[i] = ws[i] + w_grads_scale[i]
        self.w.set_weights(ws)
        print(w_grads_scale)
        print(self.w.get_weights())

        # Then, update the atomic structures.
        for i in range(len(self.hs)):
            sel = IMP.atom.Selection(self.hs[i])
            pids = sel.get_selected_particle_indexes()
            for pid in pids:
                at = IMP.atom.Atom(self.m, pid)
                at.set_occupancy(self.w.get_weight(i))





