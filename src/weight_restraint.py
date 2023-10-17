import numpy as np
import IMP
import IMP.core
import IMP.atom
import mmtbx

import xray_struct

class WeightRestraint(IMP.Restraint):
    def __init__(
            self,
            m,
            hs,
            pids,
            w,
            f_obs,
            flags,
            scale
    ):
        IMP.Restraint.__init__(self, m, "weight")
        self.m = m
        self.hs = hs
        self.pids = pids
        self.w = w
        self.n_state = len(self.hs)
        self.f_obs = f_obs
        self.flags = flags
        self.scale = scale

        self.best_f = np.infty
        self.best_occs = [self.w.get_weight(i) for i in range(self.n_state)]

    def get_scale(
            self
    ):
        return self.scale

    def get_best_occs(
            self
    ):
        return self.best_occs

    def get_w(
            self
    ):
        return self.w

    def set_scale(
            self,
            scale
    ):
        self.scale = scale

    def reset(
            self
    ):
        self.best_f = np.infty
        self.best_occs = [self.w.get_weight(i) for i in range(self.n_state)]

    def do_add_score_and_derivatives(
            self,
            sa
    ):
        # Update all occupancies based on weight of first pid.
        # n_atoms = len(self.pids)/self.n_state
        for i in range(self.n_state):
            occ = self.w.get_weight(i)
            pids = IMP.atom.Selection(self.hs[i]).get_selected_particle_indexes()
            for pid in pids:
                IMP.atom.Atom(self.m, pid).set_occupancy(occ)

            print(occ)

        xray_structure = xray_struct.get_xray_structure(
            m=self.m,
            pids=self.pids,
            crystal_symmetry=self.f_obs.crystal_symmetry()
        )

        xray_structure.scatterers().flags_set_grads(
            state=False
        )
        xray_structure.scatterers().flags_set_grad_site(
            iselection=xray_structure.all_selection().iselection()
        )
        xray_structure.scatterers().flags_set_grad_occupancy(
            iselection=xray_structure.all_selection().iselection()
        )

        f_model_manager = mmtbx.f_model.manager(
            xray_structure=xray_structure,
            f_obs=self.f_obs,
            r_free_flags=self.flags,
            target_name="ml"
        )
        f_model_manager.update_all_scales()
        fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)

        fmodels.create_target_functors()
        fmodels.prepare_target_functors_for_minimization()

        fmodels.update_xray_structure(update_f_calc = True)
        fmodels_target_and_gradients = fmodels.target_and_gradients(
            weights=None,
            compute_gradients=True,
            occupancy=True
        )

        f = fmodels_target_and_gradients.target()
        if f < self.best_f:
            self.best_f = f
            self.best_occs = [self.w.get_weight(i) for i in range(self.n_state)]

        g =  fmodels_target_and_gradients.gradients()
        occ_gs = list()

        n_atoms = len(IMP.atom.Selection(self.hs[0]).get_selected_particle_indexes())
        for i in range(self.n_state):
            occ_gs.append(g[i*n_atoms:(i+1)*n_atoms])

        avg_gs = list()
        # Get the average.
        for i in range(self.n_state):
            # avg_g = sum(occ_gs[i])/n_atoms
            avg_g = sum(occ_gs[i])
            avg_gs.append(avg_g)

        # Normalize against the average gradient.
        state_gs = list()
        for i in range(self.n_state):
            state_gs.append(avg_gs[i]-avg_gs[-1])

        print(state_gs)
        sa.add_score(f)
        if sa.get_derivative_accumulator():
            for i in range(self.n_state):
                self.w.add_to_weight_derivative(i, self.get_scale()*state_gs[i], sa.get_derivative_accumulator())

        print(self.w.get_weights_derivatives())

    def do_get_inputs(self):
        w_pid = self.w.get_particle_index()
        return [self.get_model().get_particle(w_pid)]
