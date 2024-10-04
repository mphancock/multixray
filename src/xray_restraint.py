import IMP
import IMP.core
import IMP.algebra

import miller_ops
import cctbx_score
import update_weights_optimizer_state
from derivatives import get_df_mag_ratio

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
            update_freq,
            r_charmm=None,
            ref_com=None
    ):
        IMP.Restraint.__init__(self, msmc_m.get_m(), "XrayRestraint%1%")
        self.msmc_m = msmc_m
        self.hs = msmc_m.get_hs()
        self.n_state = len(self.hs)
        self.pids = msmc_m.get_pids()
        self.cond = cond
        self.update_freq = update_freq
        self.n_evals = 0

        self.f_obs = f_obs
        self.free_flags = free_flags
        self.update_scale = update_scale
        self.update_k1 = update_k1

        self.d_min = 0
        self.set_d_min(d_min=self.d_min)

        self.target = "ml"
        self.w_xray = w_xray
        self.r_charmm = r_charmm

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
                msmc_m=self.msmc_m,
                cond=self.cond,
                f_obs=self.f_obs_filt,
                r_free_flags=self.flags_filt,
                target=self.target,
                update_scale=self.update_scale,
                update_k1=self.update_k1,
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

            self.score = score
            self.r_free = r_free
            self.r_work = r_work
            self.r_all = r_all

            if sa.get_derivative_accumulator():

                ## create dictionary for the derivatives
                ## derivatives are stored in the order of atoms/scatterers
                atoms = self.msmc_m.get_atoms()
                for i in range(len(atoms)):
                    atom = atoms[i]
                    atom_grads = grads_site[i]
                    atom_grads = IMP.algebra.Vector3D(atom_grads[0], atom_grads[1], atom_grads[2])

                    ## find the corresponding pid
                    pid = self.msmc_m.get_pid(atom)

                    # pid = self.pids[i]
                    d = IMP.core.XYZR(self.get_model(), pid)
                    dff_dx_dict[pid] = d.get_derivatives()
                    dff_avg_mag = dff_avg_mag + d.get_derivatives().get_magnitude()

                    # Store the derivative.
                    # df_dx_vec_3d = IMP.algebra.Vector3D(grads_site[i][0], grads_site[i][1], grads_site[i][2])
                    # self.df_dxs[pid] = df_dx_vec_3d

                    self.df_dxs[pid] = atom_grads

                    dxray_avg_mag = dxray_avg_mag + atom_grads.get_magnitude()


                ## this needs to be called after df dict is built
                if self.r_charmm:

                #     shadow_restraint = None
                #     for r in self.rset_charmm.get_restraints():
                #         if r.get_name() == "CHARMMShadowRestraint":
                #             shadow_restraint = r
                #             break

                    ## use the shadow derivative
                    self.df_mag_ratio = get_df_mag_ratio(
                        m=self.get_model(),
                        pids=self.pids,
                        r1=self.r_charmm,
                        r2=self
                    )
                    self.w_xray = self.df_mag_ratio
                else:
                    self.df_mag_ratio = 1

                # print(self.w_xray)

                for pid in self.pids:
                    d = IMP.core.XYZR(self.get_model(), pid)
                    dxray_dx_scaled = self.w_xray * self.df_dxs[pid]

                    d.add_to_derivatives(dxray_dx_scaled, sa.get_derivative_accumulator())

        sa.add_score(self.score)
        self.n_evals += 1
        # Store the score.


    def do_get_inputs(self):
        return [self.get_model().get_particle(pid) for pid in self.pids]

