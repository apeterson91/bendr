ITER <- 10
WARM <- 5
r <- school_data$distances
capture.output(fit <- nd_nhpp_fixed(distances_col = "distances",
                                    id_col = "school_id",
                                    data = school_data,
                                    L = 2, K = 2,
                                    iter_max = 30,
                                    warm_up = 5,
                                    thin = 1,
                                    seed = 34143))

test_that("nd_nhpp function errors if given incorrect inputs", {
    expect_error(nd_nhpp(distances_col =  "distances",
                         id_col = "school_id",
                         data = school_data,
                         iter_max = WARM, warm_up = ITER))
    expect_error(nd_nhpp_fixed(distances_col =  "distances",
                               id_col = "school_id",
                               data = school_data,
                               iter_max = WARM,
                               warm_up = ITER))
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
