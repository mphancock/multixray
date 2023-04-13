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
            pids,
            uc_dim,
            sg_symbol,
            f_obs_file,
            d_min,
            d_max,
            scale,
            target,
            w_xray,
            dynamic_w
    ):
        IMP.Restraint.__init__(self, m, "xray")
        self.pids = pids
        self.uc_dim = uc_dim
        self.sg_symbol = sg_symbol
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

        self.df_dxs = dict()
        for pid in pids:
            self.df_dxs[pid] = IMP.algebra.Vector3D(0,0,0)

        self.score = 0

        self.xray_structure = xray_struct.get_xray_structure(
            m=self.get_model(),
            uc_dim=self.uc_dim,
            sg_symbol=self.sg_symbol
        )

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

    def get_df_dict(self):
        return self.df_dxs

    def get_f(self):
        return self.score

    def do_add_score_and_derivatives(self, sa):
        # get_df_mag_ratio returns the ratio of the magnitudes of r1 to r2.
        # print("EVAL")

        # score, df_dx = cctbx_scores.get_score(
        #     m=self.get_model(),
        #     uc_dim=self.uc_dim,
        #     sg_symbol=self.sg_symbol,
        #     cif_file=self.f_obs_file,
        #     res=self.d_min,
        #     target_name=self.target
        # )

        score, df_dx = cctbx_score.get_score(
            m=self.get_model(),
            uc_dim=self.uc_dim,
            sg_symbol=self.sg_symbol,
            f_obs=self.f_obs,
            r_free_flags=self.flags,
            target=self.target
        )

        dff_dx_dict = dict()
        dff_avg_mag = 0
        dxray_avg_mag = 0
        for i in range(len(self.pids)):
            pid = self.pids[i]
            d = IMP.core.XYZR(self.get_model(), pid)
            dff_dx_dict[pid] = d.get_derivatives()
            dff_avg_mag = dff_avg_mag + d.get_derivatives().get_magnitude()

            if sa.get_derivative_accumulator():
                # df_dx_vec_3d = IMP.algebra.Vector3D(df_dx[i][0], df_dx[i][1], df_dx[i][2])
                df_dx_vec_3d = IMP.algebra.Vector3D(df_dx[i*3], df_dx[i*3+1], df_dx[i*3+2])

                # d.add_to_derivatives(df_dx_vec_3d, sa.get_derivative_accumulator())

                # Store the derivative.
                self.df_dxs[pid] = df_dx_vec_3d

                dxray_avg_mag = dxray_avg_mag + df_dx_vec_3d.get_magnitude()

        if sa.get_derivative_accumulator():
            dff_avg_mag = dff_avg_mag / len(self.pids)
            dxray_avg_mag = dxray_avg_mag / len(self.pids)

            # print(dff_avg_mag, dxray_avg_mag)
            if self.dynamic_w:
                df_mag_ratio = dff_avg_mag / dxray_avg_mag
                # print(df_mag_ratio)

            for pid in self.pids:
                w_xray = self.w_xray
                if self.dynamic_w:
                    w_xray = w_xray * df_mag_ratio

                d = IMP.core.XYZR(self.get_model(), pid)

                dxray_dx_scaled = w_xray * self.df_dxs[pid]
                d.add_to_derivatives(dxray_dx_scaled, sa.get_derivative_accumulator())

        sa.add_score(score)
        # Store the score.
        self.score = score

    def do_get_inputs(self):
        return [self.get_model().get_particle(pid) for pid in self.pids]

