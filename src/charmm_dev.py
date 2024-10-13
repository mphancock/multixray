
# class LennardJonesPairScoreWrapper(IMP.atom.LennardJonesPairScore):
#     def __init__(self, pids, sf):
#         IMP.atom.LennardJonesPairScore.__init__(self, sf)
#         self.pids = pids
#         self.df_xyz = dict()

#     def evaluate(derivs):
#         print("EVALUATE")


# class PairsRestraintWrapper(IMP.container.PairsRestraint):
#     def __init__(self, pids, score, nbl, name):
#         print("SETTING UP")
#         IMP.container.PairsRestraint.__init__(self, score, nbl, "TEST")

#         self.pids = pids
#         self.df_xyz = dict()

#     def evaluate(self, deriv):
#         print("EVALUATE")
#         return 0

#     def unprotected_evaluate(self, accum):
#         print("UNPROTECTED EVALUATE")

#     def add_score_and_derivatives(self, accum):
#         print("ADDING SCORE AND DERIVATIVES")

#     def after_evaluate(self, accum):
#         print("AFTER evaluate")
#         for pid in self.msmc_m.get_pids():
#             d = IMP.core.XYZR(self.msmc_m.get_m(), pid)
#             self.df_xyz[pid] = d.get_derivatives()

#     def do_get_inputs(self):
#         return [self.msmc_m.get_m().get_particle(pid) for pid in self.msmc_m.get_pids()]


## we can't create a wrapper around restraints (cant call evaluate or do_score_add_derviatives) so need to create a fake restraint to track the derivatives
class CHARMMRestraint(IMP.container.PairsRestraint):
    def __init__(
            self,
            msmc_m
    ):
        IMP.Restraint.__init__(self, msmc_m.get_m(), "CHARMMRestraint%1%")
        self.msmc_m = msmc_m
        self.hs = msmc_m.get_hs()
        self.n_state = len(self.hs)
        self.pids = msmc_m.get_pids()

        # Gradients and scores
        self.df_dxs = dict()
        for pid in self.pids:
            self.df_dxs[pid] = IMP.algebra.Vector3D(0,0,0)
        self.score = 0

        # setup the charmm restraints
        self.charmm_rs = list()
        for state in range(self.n_state):
            charmm_rs, nbl = charmm_restraints(
                m=self.get_model(),
                h=self.hs[state]
            )
            self.charmm_rs.extend(charmm_rs)
            self.nbl = nbl

        ## TMP IF REMOVING ATOMS
        # self.pids = IMP.atom.Selection(self.hs[0]).get_selected_particle_indexes()


    def get_df_dict(self):
        return self.df_dxs

    def get_f(self):
        return self.score

    def do_add_score_and_derivatives(self, sa):
        # self.msmc_m.get_m().do_update()
        # self.nbl.update()
        # self.nbl.do_score_state_before_evaluate()

        import time
        t0 = time.time()

        test_pid = self.pids[0]
        xyz = IMP.core.XYZR(self.get_model(), test_pid)

        self.score = 0
        for r in self.charmm_rs:
            r.add_score_and_derivatives(sa)

            print(r.get_name(), r.get_last_score())
            print("derivs", xyz.get_derivatives())

            self.score += r.get_last_score()

        print("derivs", xyz.get_derivatives())

        # create dictionary for the derivatives
        # derivatives are stored in the order of atoms/scatterers
        for pid in self.pids:
            d = IMP.core.XYZR(self.get_model(), pid)
            self.df_dxs[pid] = d.get_derivatives()

        print("NBL length", self.nbl.get_number_of_full_rebuilds())

        print("CHARMMRestraint time: ", time.time()-t0)

    def do_get_inputs(self):
        return [self.get_model().get_particle(pid) for pid in self.pids]
