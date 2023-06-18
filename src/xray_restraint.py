import IMP
import IMP.core
import IMP.algebra

import xray_struct
import miller_ops
import cctbx_score
import derivatives


class XtalRestraint(IMP.Restraint):
    def __init__(
            self,
            m,
            n_state,
            pids,
            f_obs_file,
            d_min,
            d_max,
            scale,
            target,
            w_xray,
            dynamic_w
    ):
        IMP.Restraint.__init__(self, m, "xray")
        self.n_state = n_state
        self.pids = pids
        self.f_obs_file = f_obs_file

        # Set f_obs.
        f_obs_array = miller_ops.get_miller_array(
            f_obs_file=f_obs_file,
            label="_refln.F_meas_au"
        )
        f_obs_array = miller_ops.clean_miller_array(f_obs_array)

        # Set flags.
        status_array = miller_ops.get_miller_array(
            f_obs_file=f_obs_file,
            label="_refln.status"
        )
        flags_array = status_array.customized_copy(data=status_array.data()=="f")
        f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)
        self.f_obs = f_obs_array
        self.flags = flags_array

        self.d_max = None
        self.d_min = None
        self.set_d_min(d_min=d_min)

        self.scale = scale
        self.target = target
        self.w_xray = None
        self.dynamic_w = dynamic_w

        self.set_weight(
            w_xray=w_xray
        )

        # Gradients and scores
        self.df_dxs = dict()
        self.w_grads = [0]*n_state
        for pid in pids:
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
            d_max=self.d_max,
            d_min=self.d_min
        )

    def set_d_max(
            self,
            d_max
    ):
        self.d_max = d_max
        self.f_obs_filt = miller_ops.filter_f_obs_resolution(
            f_obs=self.f_obs,
            d_max=self.d_max,
            d_min=self.d_min
        )

    def set_weight(
            self,
            w_xray
    ):
        self.w_xray = w_xray

    def set_dynamic_w(
            self,
            dynamic_w
    ):
        self.dynamic_w = dynamic_w

    def get_weight(
            self
    ):
        return self.w_xray

    def get_dynamic_w(
            self
    ):
        return self.dynamic_w

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

    def do_add_score_and_derivatives(self, sa):
        results_dict = cctbx_score.get_score(
            m=self.get_model(),
            # uc_dim=self.uc_dim,
            # sg_symbol=self.sg_symbol,
            f_obs=self.f_obs,
            r_free_flags=self.flags,
            target=self.target
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
                df_dx_vec_3d = IMP.algebra.Vector3D(grads_site[i][0], grads_site[i][1], grads_site[i][2])

                # Store the derivative.
                self.df_dxs[pid] = df_dx_vec_3d

                dxray_avg_mag = dxray_avg_mag + df_dx_vec_3d.get_magnitude()

        # Calculate and save the weights gradients.
        n_atoms = len(grads_occ) // self.n_state
        for i in range(self.n_state):
            state_grad_occs = grads_occ[i*n_atoms:(i+1)*n_atoms]
            state_grad_w = state_grad_occs[0]
            self.w_grads[i] = state_grad_w

        if sa.get_derivative_accumulator():
            dff_avg_mag = dff_avg_mag / len(self.pids)
            dxray_avg_mag = dxray_avg_mag / len(self.pids)

            if self.dynamic_w:
                df_mag_ratio = dff_avg_mag / dxray_avg_mag

            for pid in self.pids:
                w_xray = self.w_xray
                if self.dynamic_w:
                    w_xray = w_xray * df_mag_ratio

                d = IMP.core.XYZR(self.get_model(), pid)

                dxray_dx_scaled = w_xray * self.df_dxs[pid]

                # if pid == self.pids[0]:
                #     print(dxray_dx_scaled.get_magnitude())
                d.add_to_derivatives(dxray_dx_scaled, sa.get_derivative_accumulator())

        sa.add_score(score)
        # Store the score.
        self.score = score
        self.r_free = r_free
        self.r_work = r_work
        self.r_all = r_all

    def do_get_inputs(self):
        return [self.get_model().get_particle(pid) for pid in self.pids]

