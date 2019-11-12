ITER <- 10
WARM <- 5
r <- school_data %>% dplyr::arrange(school_id) %>%  ## sort to ensure correspondence with n_j
    dplyr::select(distances) %>% dplyr::pull()

n_j <- school_data %>% dplyr::arrange(school_id) %>%
    dplyr::group_by(school_id) %>% dplyr::count() %>%
    dplyr::ungroup() %>% dplyr::mutate(start = (cumsum(n) ) ) %>%
    dplyr::mutate(start_ = tidyr::replace_na(dplyr::lag(start),0) ) %>% dplyr::select(-start) %>%
    dplyr::rename(start=start_,go =n) %>%
    dplyr::select(start,go) %>% as.matrix()
capture.output(fit <- nd_nhpp_fixed(r = r, n_j = n_j,
                     L = 2, K = 2,
                     iter_max = 10,
                     warm_up = 5,
                     thin = 1,
                     seed = 34143))

test_that("nd_nhpp function errors if given incorrect inputs", {
    expect_error(nd_nhpp(c(-1,1,1),n_j = matrix(1,2,2),
                         iter_max = ITER, warm_up = WARM,),
                         regexp = "must be positive")
    expect_error(nd_nhpp(c(-1,1,1),n_j = matrix(1,2,2),
                         iter_max = WARM, warm_up = ITER),
                 regexp = "must be <")
    expect_error(nd_nhpp_fixed(c(-1,1,1),n_j = matrix(1,2,2),
                         iter_max = WARM, warm_up = ITER),
                 regexp = "must be <")
})

test_that("nd_nhpp plotting methods work",{
    expect_silent(plot(fit))
    expect_silent(plot_global_density(fit))
    expect_silent(plot_global_density(fit,r=r))
    expect_silent(plot_cluster_densities(fit))
    # expect_silent(plot_traceplots(fit)) errors on covr tests for some reason
    expect_silent(plot_pairs(fit))
    expect_silent(plot_pairs(fit,sort=TRUE))
    expect_silent(plot_network(fit,mode_label = FALSE))
    expect_silent(plot_network(fit,mode_label = TRUE))
})

test_that("ndp print and summary methods work",{
    expect_output(print(fit),"Cluster")
    expect_output(print(fit),"tau")
    expect_output(print(summary(fit)),"tau")
    expect_is(summary(fit),"summary.ndp")
})


test_that("assign error ndp methods work",{
    expect_silent(assign_mode(fit))
    expect_silent(get_square_error(fit))
    expect_silent(green_loss(fit))
})

test_that("components are mcmc objects",{
    expect_is(fit$alpha,"mcmc.list")
    expect_is(fit$rho,"mcmc.list")
    expect_is(fit$mu,"mcmc.list")
    expect_is(fit$tau,"mcmc.list")
})
