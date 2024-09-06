import numpy as np
import random

import IMP
import IMP.core
import IMP.atom

import weights


class UpdateWeightsOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            msmc_m,
            r_xrays,
            n_proposals,
            radius,
            write=True
    ):
        IMP.OptimizerState.__init__(self, msmc_m.get_m(), "UpdateWeightsOptimizerState%1%")
        self.msmc_m = msmc_m
        self.m = msmc_m.get_m()
        self.hs = msmc_m.get_hs()
        self.r_xrays = r_xrays
        self.n_proposals = n_proposals
        self.radius = radius
        self.write = write

        self.on = True

    def get_on(self):
        return self.on

    def set_on(
        self,
        on
    ):
        self.on = on

    def do_update(self, call):
        for i in range(len(self.r_xrays)):
            cur_occs = self.msmc_m.get_w_mat()[:,i]

            new_occs = self.get_new_occs(
                cur_occs=cur_occs,
                cond=i
            )

            self.msmc_m.set_occs_for_condition_i(new_occs, i)


    def get_new_occs(
        self,
        cur_occs,
        cond
    ):
        r_xray = self.r_xrays[cond]
        cur_score = r_xray.get_f()

        all_tmp_occs = [cur_occs]
        for _ in range(self.n_proposals):
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
            # print(self.get_name(), " evaluation")

            # Update the weights.
            self.msmc_m.set_occs_for_condition_i(tmp_occs, cond)

            # Don't need to compute derivatives.
            r_xray.evaluate(False)
            tmp_score = r_xray.get_f()

            if tmp_score < best_score:
                best_occs = tmp_occs
                best_score = tmp_score

            if self.write:
                print(tmp_occs, tmp_score, best_occs, best_score)

        if self.write:
            print(cur_score, cur_occs)
            print(best_score, best_occs)

        return best_occs


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



