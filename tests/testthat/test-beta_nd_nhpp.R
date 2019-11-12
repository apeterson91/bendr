ITER <- 10
WARM_UP <- 5
r <- school_data %>% dplyr::arrange(school_id) %>%  ## sort to ensure correspondence with n_j
    dplyr::select(distances) %>% dplyr::pull()

n_j <- school_data %>% dplyr::arrange(school_id) %>%
    dplyr::group_by(school_id) %>% dplyr::count() %>%
    dplyr::ungroup() %>% dplyr::mutate(start = (cumsum(n) ) ) %>%
    dplyr::mutate(start_ = tidyr::replace_na(dplyr::lag(start),0) ) %>% dplyr::select(-start) %>%
    dplyr::rename(start=start_,go =n) %>%
    dplyr::select(start,go) %>% as.matrix()

capture_output(fit <- beta_nd_nhpp(r = r, n_j = n_j,
                                   tau_sd = 1, mu_sd = 1,
                                   L = 2, K = 2,
                                   iter_max = 30,
                                   warm_up = 5,
                                   thin = 1,
                                   seed = 34143))

test_that("beta_nd_nhpp errors appropriately", {
    expect_error(beta_nd_nhpp(c(-1,1,1),n_j = matrix(1,2,2),
                         iter_max = ITER, warm_up = WARM_UP),
                 regexp = "must be positive")
    expect_error(beta_nd_nhpp(c(-1,1,0),n_j = matrix(1,2,2),
                              iter_max = ITER, warm_up = WARM_UP),
                 regexp = "must be positive")
    expect_error(beta_nd_nhpp(c(-1,1,1),n_j = matrix(1,2,2),
                              iter_max = WARM_UP, warm_up = ITER))
})

test_that("beta nd_nhpp plotting methods work",{
    expect_silent(plot(fit))
    expect_silent(plot_global_density(fit))
    expect_silent(plot_global_density(fit,r=r))
    expect_silent(plot_cluster_densities(fit))
    # expect_is(plot_traceplots(fit),"ggplot") errors on covr tests but not standard tests for some reason
    expect_silent(plot_pairs(fit))
    expect_silent(plot_pairs(fit,sort=TRUE))
    expect_silent(plot_network(fit,mode_label = FALSE))
    expect_silent(plot_network(fit,mode_label = TRUE))
})

test_that("bndp print and summary methods work",{
    expect_output(print(fit),"Cluster")
    expect_output(print(fit),"tau")
    capture_output(print(summary(fit)),"Geweke")
    expect_is(summary(fit),"summary.bndp")
})


test_that("bndp assign error ndp methods work",{
    expect_silent(assign_mode(fit))
    expect_silent(get_square_error(fit))
    expect_silent(green_loss(fit))
})

test_that("bndp components are mcmc objects",{
    expect_is(fit$alpha,"mcmc.list")
    expect_is(fit$rho,"mcmc.list")
    expect_is(fit$mu,"mcmc.list")
    expect_is(fit$tau,"mcmc.list")
})
