import numpy as np
import random

import IMP
import IMP.core
import IMP.atom

import weights


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
        cur_occs = self.w.get_weights()
        cur_score = self.r_xray.get_f()

        all_tmp_occs = [cur_occs]
        for i in range(self.n_proposals):
            tmp_occs = weights.get_weights(
                floor=.05,
                n_state=len(cur_occs),
                occs_cur=cur_occs,
                sigma=.05
            )
            all_tmp_occs.append(tmp_occs)

        # Check the score for all proposed weights.
        best_occs = cur_occs
        best_score = cur_score
        for tmp_occs in all_tmp_occs:
            print(self.get_name(), " evaluation")
            self.w.set_weights(tmp_occs)

            # Don't need to compute derivatives.
            self.r_xray.evaluate(False)
            tmp_score = self.r_xray.get_f()

            if tmp_score < best_score:
                best_occs = tmp_occs
                best_score = tmp_score

            print(tmp_occs, tmp_score, best_occs, best_score)

        self.w.set_weights(best_occs)


# class OptimizeWeightsOptimizerState(IMP.OptimizerState):
#     def __init__(
#         self,
#         m,
#         wrs,
#         w,
#         n_state,
#         step_tracker
#     ):
#         IMP.OptimizerState.__init__(self, m, "OptimizeWeightsOptimizerState%1%")
#         self.m = m
#         self.wrs = wrs
#         self.n_state = n_state
#         self.step_tracker = step_tracker

#         self.w = w
#         self.on = True

#     def get_on(self):
#         return self.on

#     def set_on(
#         self,
#         on
#     ):
#         self.on = on

#     def do_update(self, call):
#         if self.on and self.step_tracker.get_step() > 0:
#         # if self.on and self.r_xray.get_r_work() < .2:
#             for wr in self.wrs:
#                 wr.reset()

#             sf = IMP.core.RestraintsScoringFunction(self.wrs)
#             # sf.evaluate(True)

#             cg = IMP.core.ConjugateGradients(self.m)
#             cg.set_scoring_function(sf)
#             cg.optimize(10)

#             best_occs = self.wr.get_best_occs()
#             self.w.set_weights(best_occs)
#         # else:
#         #     self.w.set_weights([1/self.n_state]*self.n_state)



