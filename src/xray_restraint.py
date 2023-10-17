import time
import random
import numpy as np

import IMP
import IMP.core
import IMP.algebra

from cctbx.array_family import flex

import xray_struct
import miller_ops
import cctbx_score
import derivatives


"""
pids may only be a subset of the total pids (eg, only the pids in the main chain).
"""
class XtalRestraint(IMP.Restraint):
    def __init__(
            self,
            hs,
            pids,
            n_state,
            f_obs_file,
            scale,
            target,
            w_xray,
            dynamic_w,
            u_aniso_file,
            dropout=False,
            d_min=None,
            rand_noise=False
    ):
        IMP.Restraint.__init__(self, hs[0].get_model(), "XrayRestraint%1%")
        self.hs = hs
        self.n_state = n_state
        self.pids = pids
        self.f_obs_file = f_obs_file
        self.u_aniso_file = u_aniso_file

        # Set f_obs.
        f_obs_array = miller_ops.get_miller_array(
            f_obs_file=f_obs_file,
            label="_refln.F_meas_au"
        )
        f_obs_array = miller_ops.clean_miller_array(f_obs_array)

        rand_flags = False
        if rand_flags:
            flags_array = f_obs_array.generate_r_free_flags(
                fraction=0.2,
                max_free=len(f_obs_array.data())
            )
        else:
            # Set flags from file.
            status_array = miller_ops.get_miller_array(
                f_obs_file=f_obs_file,
                label="_refln.status"
            )
            flags_array = status_array.customized_copy(data=status_array.data()=="f")
        f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

        if d_min:
            d_filter = list()
            ds = f_obs_array.d_spacings()
            for i in range(len(ds.data())):
                # Keep all free reflections.

                if flags_array.data()[i]:
                    d_filter.append(True)
                elif ds.data()[i] > d_min:
                    d_filter.append(True)
                else:
                    d_filter.append(False)

            f_obs_array = f_obs_array.select(flex.bool(d_filter))

        if rand_noise:
            # Add noise to the work reflections.
            for i in range(len(f_obs_array.data())):
                f_ob = f_obs_array.data()[i]
                f_ob_err = np.random.normal(loc=f_ob, scale=f_ob*.05)
                f_obs_array.data()[i] = f_ob_err

        # Dropout some work reflections.
        if dropout:
            rand_select = list()

            for i in range(len(f_obs_array.data())):
                if not flags_array.data()[i]:
                    rand_select.append(random.choices([True, False], weights=[0.8, 0.2])[0])
                else:
                    rand_select.append(True)

            f_obs_array = f_obs_array.select(flex.bool(rand_select))

        # flags_array = status_array.customized_copy(data=status_array.data()=="f")

        f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)
        print("N_OBS: ", len(f_obs_array.data()), len(flags_array.data()))

        self.f_obs = f_obs_array
        self.flags = flags_array
        self.f_obs_filt = f_obs_array
        self.f_flags_filt = flags_array

        self.d_max = None
        self.d_min = None
        self.set_d_min(d_min=d_min)

        self.scale = scale
        self.target = target
        self.w_xray = None
        self.dynamic_w = dynamic_w
        self.df_mag_ratio = None

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
        self.flags_filt = miller_ops.filter_f_obs_resolution(
            f_obs=self.flags,
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

    def get_weight(self):
        return self.w_xray

    def get_dynamic_w(self):
        return self.dynamic_w

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

    def do_add_score_and_derivatives(self, sa):
        at = IMP.atom.Atom(self.get_model(), self.pids[0])

        # Get the derivatives.
        results_dict = cctbx_score.get_score(
            hs=self.hs,
            pids=self.pids,
            f_obs=self.f_obs_filt,
            r_free_flags=self.flags_filt,
            target=self.target,
            u_aniso_file=self.u_aniso_file
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
            # state_grad_w = state_grad_occs[0]
            state_grad_w = sum(state_grad_occs)
            self.w_grads[i] = state_grad_w

        if sa.get_derivative_accumulator():
            dff_avg_mag = dff_avg_mag / len(self.pids)
            dxray_avg_mag = dxray_avg_mag / len(self.pids)

            if self.dynamic_w:
                df_mag_ratio = dff_avg_mag / dxray_avg_mag
                self.df_mag_ratio = df_mag_ratio

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

        # print(self.score)

    def do_get_inputs(self):
        return [self.get_model().get_particle(pid) for pid in self.pids]

