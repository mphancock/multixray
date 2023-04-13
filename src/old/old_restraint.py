
# def do_add_score_and_derivatives(self, sa):
#     self.xray_structure = get_xray_structure(
#         m=self.get_model(),
#         uc_dim=self.uc_dim,
#         sg_symbol=self.sg_symbol
#     )
#
#     self.xray_structure.scatterers().flags_set_grads(
#         state=False
#     )
#     self.xray_structure.scatterers().flags_set_grad_site(
#         iselection=self.xray_structure.all_selection().iselection()
#     )
#
#     score, df_dx = score_and_derivatives(
#         xray_structure=self.xray_structure,
#         f_obs=self.f_obs_filt,
#         scale=self.scale,
#         target=self.target
#     )
#
#     for i in range(len(self.pids)):
#         pid = self.pids[i]
#         d = IMP.core.XYZR(self.get_model(), pid)
#
#         if sa.get_derivative_accumulator():
#             # df_dx_vec_3d = IMP.algebra.Vector3D(df_dx[i][0], df_dx[i][1], df_dx[i][2])
#             df_dx_vec_3d = IMP.algebra.Vector3D(df_dx[i*3], df_dx[i*3+1], df_dx[i*3+2])
#
#             d.add_to_derivatives(df_dx_vec_3d, sa.get_derivative_accumulator())
#
#             # Store the derivative.
#             self.df_dxs[pid] = df_dx_vec_3d
#
#     sa.add_score(score)
#     # Store the score.
#     self.score = score
