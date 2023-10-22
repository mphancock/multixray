import numpy as np
import random

import IMP
import IMP.core
import IMP.atom


def update_multi_state_model(
        hs,
        m,
        w
):
    # Update the occupancies of the model based on the weight particle attatched to first pid.
    for i in range(len(hs)):
        occ = w.get_weight(i)
        pids = IMP.atom.Selection(hs[i]).get_selected_particle_indexes()
        for pid in pids:
            IMP.atom.Atom(m, pid).set_occupancy(occ)


def get_weights_from_hs(hs):
    weights = list()
    for h in hs:
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        occ = IMP.atom.Atom(h.get_model(), pids[0]).get_occupancy()
        weights.append(occ)

    return weights


class UpdateWeightsOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            hs,
            w,
            r_xray,
            n_proposals,
            radius
    ):
        IMP.OptimizerState.__init__(self, m, "UpdateWeightsOptimizerState%1%")
        self.m = m
        self.hs = hs
        self.w = w
        self.r_xray = r_xray
        self.n_proposals = n_proposals
        self.radius = radius

        self.on = True

    def get_on(self):
        return self.on

    def set_on(
        self,
        on
    ):
        self.on = on

    def do_update(self, call):
        w_origs = self.w.get_weights()
        w_origs_score = self.r_xray.get_f()

        all_w_tmps = [w_origs]
        for i in range(self.n_proposals):
            w_tmps = np.random.normal(w_origs, self.radius)

            # if any proposed weight is 0, then randomly select a new set of weights.
            for w in w_tmps:
                if w < 1e-5:
                    w_tmps = [random.random()] * len(w_origs)
                    break

            w_tmps = [w/sum(w_tmps) for w in w_tmps]
            all_w_tmps.append(w_tmps)

        # Check the score for all proposed weights.
        best_ws = w_origs
        best_score = w_origs_score
        for w_tmps in all_w_tmps:
            print(self.get_name(), " evaluation")
            self.w.set_weights(w_tmps)

            # Don't need to compute derivatives.
            self.r_xray.evaluate(False)
            w_tmps_score = self.r_xray.get_f()

            if w_tmps_score < best_score:
                best_ws = w_tmps
                best_score = w_tmps_score

            print(w_tmps, w_tmps_score, best_ws, best_score)

        self.w.set_weights(best_ws)


class OptimizeWeightsOptimizerState(IMP.OptimizerState):
    def __init__(
        self,
        m,
        wrs,
        w,
        n_state,
        step_tracker
    ):
        IMP.OptimizerState.__init__(self, m, "OptimizeWeightsOptimizerState%1%")
        self.m = m
        self.wrs = wrs
        self.n_state = n_state
        self.step_tracker = step_tracker

        self.w = w
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
            for wr in self.wrs:
                wr.reset()

            sf = IMP.core.RestraintsScoringFunction(self.wrs)
            # sf.evaluate(True)

            cg = IMP.core.ConjugateGradients(self.m)
            cg.set_scoring_function(sf)
            cg.optimize(10)

            best_occs = self.wr.get_best_occs()
            self.w.set_weights(best_occs)
        # else:
        #     self.w.set_weights([1/self.n_state]*self.n_state)



