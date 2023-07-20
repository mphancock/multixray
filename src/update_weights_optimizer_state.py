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
        print(w_grads)

        # Normalize
        for i in range(len(w_grads)):
            w_grads[i] = w_grads[i]-w_grads[-1]

        # w_grads_scale = [w_grads[i]*1000 for i in range(len(w_grads))]
        ws = self.w.get_weights()
        for i in range(len(ws)):
            ws[i] = ws[i] + w_grads[i]
        self.w.set_weights(ws)
        print(w_grads)
        print(self.w.get_weights())

        # Then, update the atomic structures.
        for i in range(len(self.hs)):
            sel = IMP.atom.Selection(self.hs[i])
            pids = sel.get_selected_particle_indexes()
            for pid in pids:
                at = IMP.atom.Atom(self.m, pid)
                at.set_occupancy(self.w.get_weight(i))


class OptimizeWeightsOptimizerState(IMP.OptimizerState):
    def __init__(
        self,
        m,
        wr,
        r_xray,
        n_state,
        step_tracker
    ):
        IMP.OptimizerState.__init__(self, m, "OptimizeWeightsOptimizerState%1%")
        self.m = m
        self.wr = wr
        self.r_xray = r_xray
        self.n_state = n_state
        self.step_tracker = step_tracker

        self.w = self.wr.get_w()
        self.on = True

    def get_on(self):
        return self.on

    def set_on(
        self,
        on
    ):
        self.on = on

    def do_update(self, call):
        if self.on and self.step_tracker.get_step() > 0:
        # if self.on and self.r_xray.get_r_work() < .2:
            self.wr.reset()
            sf = IMP.core.RestraintsScoringFunction([self.wr])
            # sf.evaluate(True)

            cg = IMP.core.ConjugateGradients(self.m)
            cg.set_scoring_function(sf)
            cg.optimize(10)

            best_occs = self.wr.get_best_occs()
            self.w.set_weights(best_occs)
        # else:
        #     self.w.set_weights([1/self.n_state]*self.n_state)



