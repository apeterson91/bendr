ITER <- 10

capture_output(fit  <- bend(school_id ~ FFR,
                                  benvo = school_data,
                                  L = 5, K = 5,
                                  base_measure = normal_measure(),
                                  iter_max = ITER, ## To get good resolution, could possibly use more or less depending on convergence
                                  thin = 1,
                                  fix_concentration = TRUE, ## To avoid collapsing which sometimes occurs with simulated data
                                  seed = 34143))
capture_output(
    fit2 <- bend(school_id ~ FFR,
                 benvo = school_data,
                 L = 5, K = 5,
                 base_measure = beta_measure(),
                 iter_max = ITER, ## To get good resolution, could possibly use more or less depending on convergence
                 thin = 1,
                 fix_concentration = TRUE, ## To avoid collapsing which sometimes occurs with simulated data
                 seed = 34143)
)

test_that("`bend()`  errors if given incorrect inputs", {
    expect_error(bend(school_id ~ FFR,
                      benvo = school_data,
                      L = 5, K = 5,
                      iter_max = WARM, burn_in = ITER))
})

test_that("bend plotting methods work",{
    expect_silent(plot(fit))
    expect_silent(plot(fit,plotfun='global'))
    expect_silent(plot(fit,'trace'))
    expect_silent(plot(fit,'pairs'))
    expect_silent(plot(fit,'pairs',sort=TRUE))
    expect_silent(plot(fit2))
    expect_silent(plot(fit2,plotfun='global'))
    expect_silent(plot(fit2,'trace'))
    expect_silent(plot(fit2,'pairs'))
    expect_silent(plot(fit2,'pairs',sort=TRUE))
})

test_that("ndp  methods work",{
    expect_output(print(fit),"normal")
    expect_is(summary(fit),"summary.ndp")
    expect_equal(5,nsamples(fit))
    expect_output(print(fit2),"beta")
    expect_is(summary(fit2),"summary.ndp")
    expect_equal(5,nsamples(fit2))
    expect_equal(82,ncol(as.matrix(fit)))
})


test_that("assign error ndp methods work",{
    expect_silent(assign_mode(fit))
    expect_silent(get_square_error(fit))
    expect_silent(green_loss(fit))
    expect_silent(assign_mode(fit2))
    expect_silent(get_square_error(fit2))
    expect_silent(green_loss(fit2))
})

