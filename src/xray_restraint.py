import IMP
import IMP.core
import IMP.algebra

import miller_ops
import cctbx_score
import update_weights_optimizer_state


"""
pids may only be a subset of the total pids (eg, only the pids in the main chain).
"""
class XtalRestraint(IMP.Restraint):
    def __init__(
            self,
            msmc_m,
            cond,
            f_obs,
            free_flags,
            w_xray,
            update_scale,
            update_k1,
            u_aniso_file,
            update_freq,
            ref_com=None
    ):
        IMP.Restraint.__init__(self, msmc_m.get_m(), "XrayRestraint%1%")
        self.msmc_m = msmc_m
        self.hs = msmc_m.get_hs()
        self.n_state = len(self.hs)
        self.pids = msmc_m.get_all_pids()
        self.cond = cond
        self.update_freq = update_freq
        self.n_evals = 0

        self.f_obs = f_obs
        self.free_flags = free_flags
        self.update_scale = update_scale
        self.update_k1 = update_k1
        self.u_aniso_file = u_aniso_file

        self.d_min = 0
        self.set_d_min(d_min=self.d_min)

        self.target = "ml"
        self.w_xray = None
        self.df_mag_ratio = None

        self.set_weight(
            w_xray=w_xray
        )

        self.ref_com = ref_com

        # Gradients and scores
        self.df_dxs = dict()
        self.w_grads = [0]*self.n_state
        for pid in self.pids:
            self.df_dxs[pid] = IMP.algebra.Vector3D(0,0,0)
        self.score = 0
        self.r_free = 0
        self.r_work = 0
        self.r_all = 0

    def set_d_min(
            self,
            d_min
    ):
        self.d_min = d_min
        self.f_obs_filt = miller_ops.filter_f_obs_resolution(
            f_obs=self.f_obs,
            d_max=None,
            d_min=self.d_min
        )
        self.flags_filt = miller_ops.filter_f_obs_resolution(
            f_obs=self.free_flags,
            d_max=None,
            d_min=self.d_min
        )

        # if ref_hs:


    def set_weight(
            self,
            w_xray
    ):
        self.w_xray = w_xray

    def set_occs(self, occs):
        self.occs = occs

    def set_update_freq(self, update_freq):
        self.update_freq = update_freq

    def get_weight(self):
        return self.w_xray

    def get_df_mag_ratio(self):
        return self.df_mag_ratio

    def get_df_dict(self):
        return self.df_dxs

    def get_w_grads(self):
        return self.w_grads

    def get_f(self):
        return self.score

    def get_r_free(self):
        return self.r_free

    def get_r_work(self):
        return self.r_work

    def get_r_all(self):
        return self.r_all

    def get_update_freq(self):
        return self.update_freq

    def do_add_score_and_derivatives(self, sa):
        # Get the reference COM.
        if self.ref_com:
            self.delta = self.ref_com.get_coordinates() - self.msmc_m.get_com().get_coordinates()
            # print("REF COM: ", self.ref_com.get_coordinates())
        else:
            self.delta = None
        # print("DELTA: ", self.delta)

        # Get the derivatives.

        if self.n_evals % self.update_freq == 0:
            results_dict = cctbx_score.get_score(
                hs=self.msmc_m.get_hs(),
                occs=self.msmc_m.get_w_mat()[:,self.cond],
                pids=self.pids,
                f_obs=self.f_obs_filt,
                r_free_flags=self.flags_filt,
                target=self.target,
                update_scale=self.update_scale,
                update_k1=self.update_k1,
                u_aniso_file=self.u_aniso_file,
                delta=self.delta
            )

            score = results_dict["score"]
            grads_site = results_dict["grads_site"]
            grads_occ = results_dict["grads_occ"]
            r_work = results_dict["r_work"]
            r_free = results_dict["r_free"]
            r_all = results_dict["r_all"]

            dff_dx_dict = dict()
            dff_avg_mag = 0
            dxray_avg_mag = 0
            for i in range(len(self.pids)):
                pid = self.pids[i]
                d = IMP.core.XYZR(self.get_model(), pid)
                dff_dx_dict[pid] = d.get_derivatives()
                dff_avg_mag = dff_avg_mag + d.get_derivatives().get_magnitude()

                if sa.get_derivative_accumulator():
                    # Store the derivative.
                    df_dx_vec_3d = IMP.algebra.Vector3D(grads_site[i][0], grads_site[i][1], grads_site[i][2])
                    self.df_dxs[pid] = df_dx_vec_3d

                    dxray_avg_mag = dxray_avg_mag + df_dx_vec_3d.get_magnitude()

            self.score = score
            self.r_free = r_free
            self.r_work = r_work
            self.r_all = r_all

        # # Calculate and save the weights gradients.
        # n_atoms = len(grads_occ) // self.n_state
        # for i in range(self.n_state):
        #     state_grad_occs = grads_occ[i*n_atoms:(i+1)*n_atoms]
        #     # state_grad_w = state_grad_occs[0]
        #     state_grad_w = sum(state_grad_occs)
        #     self.w_grads[i] = state_grad_w

        if sa.get_derivative_accumulator():
            # dff_avg_mag = dff_avg_mag / len(self.pids)
            # dxray_avg_mag = dxray_avg_mag / len(self.pids)

            # df_mag_ratio = dff_avg_mag / dxray_avg_mag
            # self.df_mag_ratio = df_mag_ratio
            df_mag_ratio = 1

            for pid in self.pids:
                d = IMP.core.XYZR(self.get_model(), pid)
                dxray_dx_scaled = self.w_xray * df_mag_ratio * self.df_dxs[pid]

                d.add_to_derivatives(dxray_dx_scaled, sa.get_derivative_accumulator())

        sa.add_score(self.score)
        self.n_evals += 1
        # Store the score.


    def do_get_inputs(self):
        return [self.get_model().get_particle(pid) for pid in self.pids]

