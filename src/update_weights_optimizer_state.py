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
            r_xrays,
            n_proposals,
            radius
    ):
        IMP.OptimizerState.__init__(self, m, "UpdateWeightsOptimizerState%1%")
        self.m = m
        self.hs = hs
        self.w = w
        self.r_xrays = r_xrays
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
        w_origs_scores = [r_xray.get_f() for r_xray in self.r_xrays]

        all_w_tmps = [w_origs]
        all_w_tmps_scores = [w_origs_scores]

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
        for w_tmps in all_w_tmps:
            self.w.set_weights(w_tmps)
            update_multi_state_model(self.hs, self.m, self.w)

            # Check satisfaction of all x_rays.
            w_tmps_scores = list()
            for r_xray in self.r_xrays:
                # Don't need to compute derivatives.
                r_xray.evaluate(False)
                f_tmp = r_xray.get_f()
                w_tmps_scores.append(f_tmp)

            all_w_tmps_scores.append(w_tmps_scores)

        best_ws = all_w_tmps[0]
        best_scores = all_w_tmps_scores[0]
        # Find the best set of weights:
        for i in range(len(all_w_tmps)):
            ws = all_w_tmps[i]
            ws_scores = all_w_tmps_scores[i]

            best = True
            for j in range(len(ws_scores)):
                if ws_scores[j] > best_scores[j]:
                    best = False
                    break

            if best:
                best_ws = ws
                best_scores = ws_scores

            print(ws, ws_scores, best)


        # print("best: ", best_ws, best_scores)
        self.w.set_weights(best_ws)
        # print(self.w.get_weights())
        update_multi_state_model(self.hs, self.m, self.w)


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



