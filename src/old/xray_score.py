import mmtbx.f_model

# def likelihood_and_derivatives(
#         xray_structure,
#         f_obs,
#         scale=True
# ):
#     f_model_manager = mmtbx.f_model.manager(
#         xray_structure=xray_structure,
#         f_obs=f_obs,
#         target_name="ml"
#     )
#
#     if scale:
#         f_model_manager.update_all_scales(
#             remove_outliers=False
#         )
#
#     functor = f_model_manager.target_functor()
#     result = functor(compute_gradients=True)
#     d_target_d_site = result.d_target_d_site_cart()
#     log_likelihood = result.target_work()
#
#     # f_model_manager.target_w()
#
#     return log_likelihood, d_target_d_site

def score_and_derivatives(
        xray_structure,
        f_obs,
        scale,
        target
):
    # f_model_manager = mmtbx.f_model.manager(
    #     xray_structure=xray_structure,
    #     f_obs=f_obs,
    #     target_name=target
    # )
    #
    # if scale:
    #     f_model_manager.update_all_scales(
    #         remove_outliers=False
    #     )
    #
    # functor = f_model_manager.target_functor()
    # result = functor(compute_gradients=True)
    # d_target_d_site = result.d_target_d_site_cart()
    # score = result.target_work()

    # f_model_manager.target_w()

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        target_name=target
    )

    fmodels = mmtbx.fmodels(fmodel_xray=f_model_manager)
    fmodels.update_xray_structure(
        xray_structure=xray_structure,
        update_f_calc=True)

    fmodels_target_and_gradients = fmodels.target_and_gradients(
        compute_gradients=True)
    target = fmodels_target_and_gradients.target()
    d_target_d_site = fmodels_target_and_gradients.gradients()

    return target, d_target_d_site
